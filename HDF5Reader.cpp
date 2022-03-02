// Open and read a file from the Acqiris system (from a single instrument)
// The H5 API is used - to compile use h5c++ instead of g++
// In order to also include root libraries:
// h5c++ `root-config --cflags` -g -O2 -o HDF55Reader HDF55Reader.cpp `root-config --libs
// Update - use the following to compile instead
// g++ `root-config --cflags` -g -O2 -o HDF5Reader HDF5Reader.cpp `root-config --libs` -I /usr/include/hdf5/serial/ -lhdf5_serial -lhdf5_cpp
// The usage is

/***************************************************/
/******************** Preamble *********************/
/***************************************************/

// c++ libraries
#include <iostream>
#include <string>
#include <vector>
#include "unistd.h"

// HDF5 libraries
#include "H5Cpp.h"

// Root libraries
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"
#include "TArray.h"
#include "TApplication.h"
#include "TText.h"

// Namespaces
using namespace H5;
using namespace std;


//***************************************************//
//******************* Functions *********************//
//***************************************************//


//*** Check if a pointer is allocated memory correctly ***//

void ptrcheck(void *ptr, string name){
  if (ptr == nullptr){
    cout << "Memory allocation failure for variable: " << name << endl;    
    exit(1);
  }
}


//*** Create calibrated and smoothed histogram from the raw array ***//

void makehist(short int *mydata, TH1F *wf_smooth, int samples, double baseline, double gain, int nsmooth){
  
  // If the smoothing factor is zero (no smoothing) just need to account for baseline and gain 
  if (nsmooth == 0){
    for (int i = 0; i < samples; i++){
      wf_smooth->SetBinContent(i+1, (mydata[i]-baseline)*gain);
    }
  }
  
  // Otherwise calculate the new value point by point
  else{
    // The histogram should be filled from bin 1, not bin 0
    // The first and last nsmooth bins should be set to zero
    for (int i = 0; i < nsmooth; i++){
      wf_smooth->SetBinContent(i+1,0);
    }
    for (int i = samples-nsmooth; i<samples; i++){
      wf_smooth->SetBinContent(i+1,0);
    }
    
    // Initial sum
    double sum = 0;
    for (int i = 0; i < (2*nsmooth)+1; i++){
      sum += mydata[i];
    }
    
    // Fill in smoothed values and calulate new sum
    for (int i = nsmooth; i<(samples-nsmooth); i++){
      double value = sum / ((2*nsmooth)+1);                  // The smoothed value
      wf_smooth->SetBinContent(i+1,(value-baseline)*gain);   // Remove the baseline and convert to mV
      sum -= mydata[i-nsmooth];
      sum += mydata[i+nsmooth+1]; 
    }
  }
}

//***************************************************//
//********************** Main ***********************//
//***************************************************//

// Expect an argument containing the h5 file to be opened
int main(int argc, char *argv[]){
  
  
  // Read arguments
  // 
  // d - Draw waveforms on event by event basis
  // o - create a root file and save TTree
  // t - define the minimum threshold for ToT calculations
  // i - the amount to increment the threshold by in ToT calculations
  // s - the smoothing factor, number of adjacent channels to add
  
  
  bool DrawWF = false;      // Draw the waveforms for each event
  bool RootOutput = false;  // Output a root file
  double MinThresh = 100;   // Minimum value for ToT calculations
  double ThreshDelta = 25;  // Increments in which the ToT thresholds are increased
  int nsmooth = 1;
  string  H5FileName = "run133.h5";
  
  
  int option;
  // Colon following an option indicates a whitespace-separated argument is expected
  while ((option = getopt (argc, argv, "h:dot:i:s:")) != -1){
    switch (option){
      case 'h':
        H5FileName = optarg;
        cout << "Input file = " << optarg << endl;
        break;
      case 'd':
        DrawWF = true;
        cout << "Waveforms will be drawn for each event" << endl;
        break;
      case 'o':
        RootOutput = true;
        cout << "An output root file will be created" << endl;
        break;
      case 't':
        MinThresh = atoi(optarg);
        cout << "The minimum ToT threshold is set to " << MinThresh << endl;
        break;
      case 'i':
        ThreshDelta = atoi(optarg);
        cout << "The ToT threshold increment is set to " << ThreshDelta << endl;
        break;
      case 's':
        nsmooth = atoi(optarg);
        cout << "The smoothing factor is set to " << nsmooth << endl;
        break;
      default:
        cout << "Unrecognised option" << endl;
        cout << "Usage: " << argv[0] << "-h inputfile.h5 [-d] [-o] [-t threshold minimum] [-i threshold increment] [-s smoothing factor]" << endl;
        cout << "-d will cause all wave forms to be drawn on an event by event basis" << endl;
        cout << "-o will cause an output root file to be created" << endl;
        exit(1);
      }
  }
  
  
  
  
  /*
  // Check arguments - expect 2
  if (argc != 3){
    cout << "Usage: " << argv[0] << " inputfile.h5 SmoothingValue" << endl;
    cout << "SmoothingValue should be an integer (try 3)" << endl;
    return 1;
  }

  // Input file name
  string  H5FileName = argv[1];
  cout << "Processing file :" << H5FileName << endl;
  
  // Smoothing factor
  int nsmooth = atoi(argv[2]);
  if (nsmooth < 1){
    cout << "The smoothing factor should be a positive integer" << endl;
    return 2;
  }*/
  

  // If any kind of graphical objects from ROOT are used, a TApplication must be created
  // Give it a fake argc and argv to stop it catching arguments being passed
  int fargc = 0;    // No arguments
  char *fargv[1]; fargv[0] = argv[0];     // The executable name
  TApplication app("app", &fargc, fargv);  

  // Variables for navigating the hdf5 file structure
  Group     group;        // A directory
  Attribute attrib;       // Metadata (value or string)
  StrType   stype;        // The data type needed for string attributes
  DataSet   dataset;      // The waveforms

  // HDF5 header information
  string  dataType;       // Type used for each sample, expect 'I8'
  string  filename;       // The file being read
  int     nAcquisitions;  // The number of triggers
  int     nChannels;      // The number of channels (total waveforms = triggers * channels)
  
  // Run parameters
  int     samples;        // The number of samples per waveform
  double  delay;          // The time to start recording relative to the trigger
  double  increment;      // The timescale

  //*****************************************************//
  //********** Open file and extract metadata ***********//
  //*****************************************************//
  
  // Open input file
  H5File* file = new (nothrow) H5File(H5FileName, H5F_ACC_RDONLY);
  if (file == nullptr){
    cout << "file" << endl;
    exit(1);    
  }
  
  // Go to header
  group  = file->openGroup("/Header");
  
  // Open various attributes from the header directory and record the values
  attrib = group.openAttribute("dataType");
  stype = attrib.getStrType();
  attrib.read(stype, dataType);  
  attrib = group.openAttribute("filename");
  stype = attrib.getStrType();
  attrib.read(stype, filename);
  attrib = group.openAttribute("nAcquisitions");
  attrib.read(PredType::STD_U32LE, &nAcquisitions);  
  attrib = group.openAttribute("nChannels");
  attrib.read(PredType::STD_U32LE, &nChannels);
  
  // Look at the run parameters
  group  = file->openGroup("/Header/RunParameters");
  attrib = group.openAttribute("delay");
  attrib.read(PredType::IEEE_F64LE, &delay);
  attrib = group.openAttribute("increment");
  attrib.read(PredType::IEEE_F64LE, &increment);
  attrib = group.openAttribute("samples");
  attrib.read(PredType::STD_U32LE, &samples);
  
  double  *gain = new (nothrow) double[nChannels];          // The gain (multiply by recorded value to get volts) - set for each channel
  ptrcheck(gain, "gain");
    
  double  *offset = new (nothrow) double[nChannels];        // The offset applied to the voltage range - set for each channel
  ptrcheck(offset, "offset");
  
  double  *baseline = new (nothrow) double[nChannels];      // Derived parameter from the gain and offset
  ptrcheck(baseline, "baseline");
  
  // Container for the waveform
  short int *mydata = new (nothrow) short int[samples];
  ptrcheck(mydata, "mydata");
  
  // Get the gain and offset for each Channel
  for (int chan=0; chan<nChannels; chan++){
    string ChanName = Form("/Header/Instrument/Channel.ch%d", chan+1);
    group  = file->openGroup(ChanName);
    attrib = group.openAttribute("gain");
    attrib.read(PredType::IEEE_F64LE, &gain[chan]); 
    attrib = group.openAttribute("offset");
    attrib.read(PredType::IEEE_F64LE, &offset[chan]);
    baseline[chan] = -offset[chan]/gain[chan]; // The expected offset of the baseline
  }

  // print metadata
  cout << "==> filename: "  << filename 
    << "\n\tnAcq: "  << nAcquisitions
    << "\n\tnCh: "     << nChannels
    << "\n\tdataType: " << dataType
    << "\n\tdelay: "   << delay
    << "\n\tinc: "     << increment
    << "\n\tsamples: " << samples
    << endl;  
  for (int chan=0; chan<nChannels; chan++){
    cout << "Channel " << chan << " gain: " << gain[chan] << endl;  
  }
  
  // Calculate the channel where the trigger should occur
  double TriggerChan = -delay/increment;
  cout << "The trigger is at channel " << TriggerChan << endl;
  
  //*****************************************************//
  //************ Set up ROOT file and TTree *************//
  //*****************************************************//  
  
  
  // Histogram array for storing waveforms
  TH1F **wf = new TH1F*[nChannels];
  for (int chan=0; chan<nChannels; chan++){
    TString histname = Form("Chan%d",chan);
    wf[chan] = new TH1F(histname,histname,samples,0,samples);
  }
  
  // Array of amplitudes
  double *Amp = new (nothrow) double [nChannels]; 
  ptrcheck(Amp, "Amp");
  
  // Arrays of leading and falling edges
  int LE[20][10];
  int FE[20][10];
  
  // If RootOutput is requested, create an output file and TTree
  // Filename is the original filename but with a .root extension
  
  TFile *outfile;
  TTree *data;
  
  if (RootOutput){
    TString RootFileName = H5FileName.erase(H5FileName.size()-2,2) + "root";
    outfile = new TFile(RootFileName,"RECREATE");
    data = new TTree("data","data");
    data->Branch("nChannels",&nChannels,"nChannels/I"); 
    data->Branch("Amp",Amp,"Amp[nChannels]/D"); 
    data->Branch("LE",LE,"LE[20][10]/I"); 
    data->Branch("FE",FE,"FE[20][10]/I");
    
    // Write some of the metadata (retrieve with e.g. samples->GetUnigueID()
    TObject nSamples;
    nSamples.SetUniqueID(samples);
    nSamples.Write("samples");
  }

  
  //*****************************************************//
  //************ Create a canvas (if wanted) ************//
  //*****************************************************//    
  
  // If the option to draw the results is not selected, this never gets drawn
  TCanvas *can;
  TPad *histPad;
  
  TPad *histTitle;
  TText *t1;
  
  if (DrawWF){
    can = new TCanvas("can","can",1200,1200);
    histPad = new TPad("Graphs","Graphs",0.01,0.05,0.95,0.95);
    histPad->Draw();
    histPad->cd();
    histPad->Divide(4,4);
    
    can->cd();
    histTitle = new TPad("Title","Title",0.1,0.96,0.9,0.99);
    histTitle->Draw();
    t1 = new TText(0.5,0.5,"Title goes here");
    t1->SetTextSize(0.8);
    t1->SetTextAlign(22);
  }
  
  
  
  //*****************************************************//
  //************ Loop over all acquisitions *************//
  //*****************************************************//   
  
  
  // Loop from Acquisition.1 to Acquisition.nAcquisitions
  for (int acqi=1; acqi<nAcquisitions+1; acqi++){  
    
    // Print status every 10 events
    if (acqi%10==0){
      cout << "Reading event " << acqi << " of " << nAcquisitions << "\r" << flush;
    }
  
    // Open acquisition directory
    group = file->openGroup(Form("/Acquisition.%d",acqi));
    
    for (int chan=0; chan<nChannels; chan++){
    
      // Open the dataset
      dataset = group.openDataSet(Form("Segment.%d.ch%d.0",acqi,chan+1));
      dataset.read(mydata, PredType::NATIVE_SHORT);
      
      // Get the smoothed histogram
      // raw array, histogram, number of samples, baseline, gain (mV), smoothing factor (number of neighbouring channels)
      makehist(mydata, wf[chan], samples, baseline[chan], gain[chan]*1000, nsmooth);
      Amp[chan] = wf[chan]->GetMaximum();


      // Get the ToT parameters - 9 different thresholds (plus one for the digital signals)
      for (int k=0; k<10; k++){
        double Threshold = MinThresh + (k * ThreshDelta);
        if (k==9) Threshold = 1500; // Set the final threshold to 1.5 V for the digital ToT signals in some runs
        LE[chan][k] = wf[chan]->FindFirstBinAbove(Threshold);
        FE[chan][k] = wf[chan]->FindLastBinAbove(Threshold);     
      }
      
    }
    
    // Fill the TTree parameters
    if (RootOutput){
      data->Fill();
    }
    
    // Draw if requested
    if (DrawWF){
      for (int chan=0; chan<16; chan++){ 
        histPad->cd(chan+1);
        wf[chan]->Draw();
      }
      
      // Update the canvas title
      histTitle->cd();
      TString TitleText = H5FileName + Form(", event %d (smoothing %d)",acqi,nsmooth);
      t1->SetText(0.5,0.5,TitleText);
      t1->Draw();
    
      can->WaitPrimitive();
    }
    
  }
  
  // Write the final TTree to file
  if (RootOutput){
    data->Write();
  }

  //*****************************************************//
  //********************* Clean up **********************//
  //*****************************************************//  
  
  // Close the h5 file
  cout << "\nClosing h5 file" << endl;      
  file->close();
  delete file;   
  
  delete [] mydata;
  delete [] gain;
  delete [] offset;
  delete [] baseline;
  delete [] Amp;
  
  // Delete root objects before closing output file
  
  // Only created if drawing
  if (DrawWF){ 
    delete can;
    delete histPad;
    delete histTitle;
    delete t1;
  }
  
  // Only created if outputting a root file
  if (RootOutput) delete data;
  
  for (int chan=0; chan<nChannels; chan++){
    delete wf[chan];
  }
  delete [] wf;

  if (RootOutput){
    cout << "Closing root file" << endl;  
    outfile->Close();
  }
  delete outfile;    
    
  cout << "Done!" << endl; 
  return 0;

}

