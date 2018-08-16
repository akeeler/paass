/** \file E14060Processor.cpp
 * \cleaned-up Experiment-specific processor for the e14060 experiment at 
 * the NSCL.  Based on a processor written by S. V. Paulauskas.
 *\ author A. Keeler
 *\ date April 10, 2018
 */

#include <iostream>
#include "DammPlotIds.hpp"

#include "DetectorDriver.hpp"
#include "GeProcessor.hpp"
#include "E14060Processor.hpp"
#include "PspmtProcessor.hpp"
#include "TimingMapBuilder.hpp"
#include "VandleProcessor.hpp"

using namespace std;
using namespace dammIds::experiment;

namespace dammIds{
  namespace experiment{
    // pin1 - pin2 comparisons
    const int DD_PIN1_PIN2 = 0;
    const int DD_PIN1_PIN2_CUT = 1;


    // PID and position correction plots
    const int DD_PIN1_TOF_PIN1_I2N = 2;
    const int DD_I2POS = 3;

    const int DD_TOF_PIN1_I2N_VS_I2NS = 4;
    const int DD_TOF_PIN1_I2N_VS_I2POS1 = 5;
    const int DD_I2NS_COR_TOF_PIN1_I2N = 6;
    const int DD_I2POS1_COR_TOF_PIN1_I2N = 7;
    const int DD_PIN1_VS_I2NS_COR_TOF_PIN1_I2N = 8;
    const int DD_PIN1_VS_I2POS1_COR_TOF_PIN1_I2N = 9;

    // Isomer ID plots
    const int DD_PIN1_GE = 10;
    const int DD_GE_COR_TOF = 11;
    const int D_GE_IMPLANT = 12;
    const int D_GE_ISOMER  = 13;

    // Decay-gated Ge
    const int D_GE_DECAY = 15;
    const int D_GE_CORR_DECAY = 16;
    const int D_GE_NON_CORR_DECAY = 17;
    

    // PSPMT diagnostics
    const int D_DYNODE_MULT = 20;
    const int D_ANODE_MULT = 21;
    const int D_BAD_ANODE_MULT = 22;
    const int DD_ANODE_DYNODE = 23;


    const int D_IMPLANT_DYNODE = 24;
    const int D_DECAY_DYNODE = 25;
    const int D_VETO_DYNODE = 26;

    const int DD_IMPLANT_POS = 27;
    const int DD_DECAY_POS = 28;
    const int DD_VETO_POS = 29;

    const int DD_IMPLANT_PIXEL = 30;
    const int DD_HOTSPOT_MAP = 31;
    const int D_DECAY_CURVE = 32;

    const int DD_IMPLANT_PIXELS = 37;
    const int DD_DECAY_PIXELS = 38;

    const int DD_XA_DYNODE = 51;
    const int DD_XB_DYNODE = 52;
    const int DD_YA_DYNODE = 53;
    const int DD_YB_DYNODE = 54;

    
  }
}

E14060Processor::~E14060Processor(){
  rootfile->Write();
  rootfile->Close();
}


void E14060Processor::DeclarePlots(void){
  DeclareHistogram2D(DD_PIN1_PIN2, SB, SB, "Pin1 vs Pin2");
  DeclareHistogram2D(DD_PIN1_PIN2_CUT, SB, SB, "Pin1 vs Pin2 - cut");

  DeclareHistogram2D(DD_PIN1_TOF_PIN1_I2N, SB, SA, "Pin1 vs TOF Pin1-I2n");
  DeclareHistogram2D(DD_I2POS, SC, S2, "I2 position from I2ns and Pin1 (I2s - I2n)");
  DeclareHistogram2D(DD_TOF_PIN1_I2N_VS_I2NS, SB, SB, "TOF Pin1-I2n vs I2ns position");
  DeclareHistogram2D(DD_TOF_PIN1_I2N_VS_I2POS1, SB, SB, "TOF Pin1-I2n vs I2s -I2n position");
  DeclareHistogram2D(DD_I2NS_COR_TOF_PIN1_I2N, SB, SB, "I2ns corrected TOF vs I2ns position");
  DeclareHistogram2D(DD_I2POS1_COR_TOF_PIN1_I2N, SB, SB, "I2s - I2n corrected TOF vs I2s - I2n Position");
  DeclareHistogram2D(DD_PIN1_VS_I2NS_COR_TOF_PIN1_I2N, SB, SA, "Pin1 vs I2ns Corrected TOF");
  DeclareHistogram2D(DD_PIN1_VS_I2POS1_COR_TOF_PIN1_I2N, SB, SA, "Pin1 vs I2s -I2n corrected TOF");

  DeclareHistogram2D(DD_PIN1_GE, SE, SA, "Pin1 vs Ge");
  DeclareHistogram2D(DD_GE_COR_TOF, SB, SE, "Ge vs Corrected ToF");
  DeclareHistogram1D(D_GE_IMPLANT, SE, "Implant-gated Ge");
  DeclareHistogram1D(D_GE_ISOMER, SE, "78Zn gated Ge");

  DeclareHistogram1D(D_GE_DECAY, SE, "Decay-gated Ge");
  DeclareHistogram1D(D_GE_CORR_DECAY, SE, "correlated decay-gated Ge");
  DeclareHistogram1D(D_GE_NON_CORR_DECAY, SE, "non-correlated Ge");

  DeclareHistogram1D(D_DYNODE_MULT, S5, "Multiplicity of implant dynode and decay dynode (+10)");
  DeclareHistogram1D(D_ANODE_MULT, S5, "Multiplicity of implant anodes and decay anodes (+10)");
  DeclareHistogram1D(D_BAD_ANODE_MULT, S5, "Multiplicity of non-position anodes, implants and decays (+10)");
  DeclareHistogram2D(DD_ANODE_DYNODE, SE, SE, "anode E sum vs dynode E");

  DeclareHistogram1D(D_IMPLANT_DYNODE, SD, "Implant-gated Low-gain dynode");
  DeclareHistogram1D(D_DECAY_DYNODE, SE, "Decay-gated Hi-gain dynode");
  DeclareHistogram1D(D_VETO_DYNODE, SD, "Veto-gated Hi-gain dynode");

  DeclareHistogram2D(DD_IMPLANT_POS, SA, SA, "Positions of implants");
  DeclareHistogram2D(DD_DECAY_POS, SA, SA, "Positions of decays");
  DeclareHistogram2D(DD_VETO_POS, SA, SA, "Positions of vetos");

  DeclareHistogram2D(DD_IMPLANT_PIXEL, SA, SA, "78Zn Implants in single pixel");
  DeclareHistogram2D(DD_HOTSPOT_MAP, S6, S6, "Implants correlated with single-pixel decays");
  DeclareHistogram1D(D_DECAY_CURVE, SB, "time difference between 1-pixel implant and decay");
  DeclareHistogram2D(DD_IMPLANT_PIXELS, SA, SA, "map of diagonal pixels for implants");
  DeclareHistogram2D(DD_DECAY_PIXELS, SA, SA, "map of diagonal pixels for decays");
  DeclareHistogram2D(DD_XA_DYNODE, SD, SD, "low-gain xa anode vs dynode");
  DeclareHistogram2D(DD_XB_DYNODE, SD, SD, "low-gain xb anode vs dynode");
  DeclareHistogram2D(DD_YA_DYNODE, SD, SD, "low-gain ya anode vs dynode");
  DeclareHistogram2D(DD_YB_DYNODE, SD, SD, "low-gain yb anode vs dynode");

  //-------------------Set ROOT Outputs----------------------

  
  /*
  roottree->Branch("Energy", &energy);
  roottree->Branch("Type", &type);
  roottree->Branch("subType", &subType);
  roottree->Branch("Tag", &tag);
  */

}

E14060Processor::E14060Processor(std::pair<double, double> &energyRange) :
  EventProcessor(OFFSET, RANGE, "E14060Processor"){
  SetAssociatedTypes();
  energyRange_ = energyRange;

  eventType = "";
  low_xa = 0;
  low_xb=0;
  low_ya=0;
  low_yb=0;
  hi_xa=0;
  hi_xb=0;
  hi_ya=0;
  hi_yb=0;
  low_dynode=0;
  hi_dynode=0;
  pos_x = 0;
  pos_y = 0;
  timestamp = 0;

  std::string fname = Globals::get()->GetOutputFileName();
  fname += ".root";
  rootfile = new TFile(fname.c_str(), "RECREATE");
  roottree = new TTree("E14060Processor","E14060Processor");
  roottree->Branch("xa_low", &low_xa);
  roottree->Branch("xa_hi", &hi_xa);
  roottree->Branch("xb_low", &low_xb);
  roottree->Branch("xb_hi", &hi_xb);
  roottree->Branch("ya_low", &low_ya);
  roottree->Branch("ya_hi", &hi_ya);
  roottree->Branch("yb_low", &low_yb);
  roottree->Branch("yb_hi", &hi_yb);

  roottree->Branch("event_type", &eventType);
  roottree->Branch("dynode_low", &low_dynode);
  roottree->Branch("dynode_hi", &hi_dynode);
  roottree->Branch("pos_x", &pos_x);
  roottree->Branch("pos_y", &pos_y);
  roottree->Branch("timestamp", &timestamp);
  //cout<<"1"<<endl;
}

void E14060Processor::SetAssociatedTypes() {
  associatedTypes.insert("vandle");
  associatedTypes.insert("hagrid");
  associatedTypes.insert("pspmt");
  associatedTypes.insert("tac");
  associatedTypes.insert("ge");
  associatedTypes.insert("si");
}

bool E14060Processor::Process(RawEvent &event) {
  if (!EventProcessor::Process(event))
    return(false);

  //---------------------Declare and Initialize values------------
      // Plotting values

  static const double plotMult = 2;
  static const double plotOffset = 1000;
  static const double pspmtScale = 500;
  static const double pspmtOffset = 500;

      //correlation values
  map<string, double> pins_and_tacs;
  pair<double, double> position = {0, 0};
  pair<int, int> pixel = {0, 0};

      //other initializations
  BarMap vbars;
  vector<ChanEvent *> geEvts;
  vector<ChanEvent *> vetoEvts;



  //-------------------Read in Event info-------------------------

  static const vector<ChanEvent *> &dynodeLow =
    event.GetSummary("generic:dynode")->GetList();
  static const vector<ChanEvent *> &dynodeHi =
    event.GetSummary("pspmt:dynode")->GetList();

  if (event.GetSummary("ge")->GetList().size() != 0)
    geEvts = event.GetSummary("ge")->GetList();
  if (event.GetSummary("generic:veto")->GetList().size() != 0)
    vetoEvts = event.GetSummary("generic:veto")->GetList();

  static const vector<ChanEvent *> &tac =
    event.GetSummary("tac", true)->GetList();
  static const vector<ChanEvent *> &pin =
    event.GetSummary("pin", true)->GetList();

  /*   
  if (event.GetSummary("vandle")->GetList().size() != 0)
    vbars = ((VandleProcessor *) DetectorDriver::get->
	     GetProcessor("VandleProcessor"))->GetBars();
  */

    //fill the pins_and_tacs map and filter into variables
  for (vector<ChanEvent *>::const_iterator it = tac.begin(); it != tac.end();
       it++){
    pins_and_tacs.insert(make_pair((*it)->GetChanID().GetSubtype(),
				   (*it)->GetCalibratedEnergy()));
  }
  for (vector<ChanEvent *>::const_iterator it = pin.begin(); it != pin.end();
       it++){
    pins_and_tacs.insert(make_pair((*it)->GetChanID().GetSubtype(),
				   (*it)->GetCalibratedEnergy()));
  }

  double pin1 = pins_and_tacs.find("de1")->second;
  double pin2 = pins_and_tacs.find("de2")->second;

  double pin1_i2n = pins_and_tacs.find("pin1_i2n")->second;
  double pin1_i2s = pins_and_tacs.find("pin1_i2s")->second;
  double pin2_i2n = pins_and_tacs.find("pin2_i2n")->second;
  double pin2_i2s = pins_and_tacs.find("pin2_i2s")->second;

  double i2ns = pins_and_tacs.find("i2n_i2s")->second;
  double i2pos1 = 1.8 * (pin1_i2s - pin1_i2n);
  double i2pos2 = 1.8 * (pin2_i2s - pin2_i2n);

  plot(DD_I2POS, i2ns, 0);
  plot(DD_I2POS, i2pos1, 1);
  plot(DD_I2POS, i2pos2, 2);

  double i2ns_cor_tof_pin1_i2n = 0;
  double i2pos1_cor_tof_pin1_i2n = 0;
  

  //-------------------- Implant/decay Logic ---------------------

  bool hasIon = !pin.empty() && pin1 > 1;
  bool hasLightIon = false;
  if(hasIon && pin1 <= 400)
    hasLightIon = true;

  bool hasVeto = false;
  bool hasHPlus = false;
  vector<ChanEvent *>::const_iterator iterator1 = vetoEvts.begin();
  while (!hasVeto &&iterator1 != vetoEvts.end()){
    if ((*iterator1)->GetCalibratedEnergy() > 400){
      hasVeto = true;
      if ((*iterator1)->GetCalibratedEnergy() < 800)
	hasHPlus = true;
    }
    iterator1++;
  }

  bool hasDynodeLow = !dynodeLow.empty() && (*dynodeLow.begin())->GetCalibratedEnergy() > 1;
  bool hasDynodeHi = !dynodeHi.empty() && (*dynodeHi.begin())->GetCalibratedEnergy() > 1;
  if(hasDynodeLow)
    low_dynode = (*dynodeLow.begin())->GetCalibratedEnergy();
  if(hasDynodeHi)
    hi_dynode = (*dynodeHi.begin())->GetCalibratedEnergy();

  bool hasImplantReject = hasLightIon || hasVeto;
  bool hasDecayReject = hasVeto || hasIon;

  bool hasImplant = hasIon && hasDynodeLow && !hasImplantReject;
  bool hasDecay = hasDynodeHi && !hasDecayReject;
  if(hasImplant)
    timestamp = (*dynodeLow.begin())->GetTimeSansCfd();

  //cout<<"2"<<endl;

  if (hasImplant)
    eventType = "implant";
  else if (hasDecay)
    eventType = "decay";
  else if (hasVeto)
    eventType = "veto";
  else
    eventType = " ";

  //cout<<"event = "<< eventType <<endl;

  bool has72Co = false;
  bool has78Zn = false;
  bool hasPID = false;
  bool has72CoDecay = false;


  //-------------------PID plotting-----------------------

  plot(DD_PIN1_PIN2, pin2, pin1);
  double delta = 0.87558748 * pin2 + 8.72618557 - pin1;
  if(pin1 != 0 && pin2 != 0 && delta < 40 && delta > -40)
    plot(DD_PIN1_PIN2_CUT, pin2, pin1);

  if(hasImplant && delta < 40 && delta > -40){

    if(pin1_i2n != 0){
      plot(DD_PIN1_TOF_PIN1_I2N, pin1_i2n, pin1);

      plot(DD_TOF_PIN1_I2N_VS_I2NS, i2ns, pin1_i2n);
      plot(DD_TOF_PIN1_I2N_VS_I2POS1, i2pos1, pin1_i2n);
      // Make position correction and plot corrected PID
      i2ns_cor_tof_pin1_i2n = CorrectTofByPos(i2ns, pin1_i2n);
      plot(DD_I2NS_COR_TOF_PIN1_I2N, i2ns, i2ns_cor_tof_pin1_i2n);
      plot(DD_PIN1_VS_I2NS_COR_TOF_PIN1_I2N, i2ns_cor_tof_pin1_i2n, pin1);
	   
      i2pos1_cor_tof_pin1_i2n = CorrectTofByPos(i2pos1, pin1_i2n);
      plot(DD_I2POS1_COR_TOF_PIN1_I2N, i2pos1, i2pos1_cor_tof_pin1_i2n);
      plot(DD_PIN1_VS_I2POS1_COR_TOF_PIN1_I2N, i2pos1_cor_tof_pin1_i2n, pin1);
      
      hasPID = true;
      
      
      if (i2pos1_cor_tof_pin1_i2n < 633 && i2pos1_cor_tof_pin1_i2n > 459){
	if (pin1 < 540 && pin1 > 504)
	  has72Co = true;
      }
      if (i2pos1_cor_tof_pin1_i2n < 1016 && i2pos1_cor_tof_pin1_i2n > 885){
	if (pin1 < 664 && pin1 > 617)
	  has78Zn = true;
      }
    }
  }

  //-------------------Calculate positions for Pspmt---------------

  double xa = 0; double xb = 0; double ya = 0; double yb = 0;
  //const int Px = 12; const int Py = 12;
  double PspmtThreshold = 50;
  bool hasPosition = false;

  vector<ChanEvent *> PSPMTanode = {};
  vector<ChanEvent *> OtherAnode = {};
  
  if(hasImplant){
    PSPMTanode = event.GetSummary("generic:anode")->GetList();
    OtherAnode = event.GetSummary("pspmt:anode")->GetList();
  }
  else if (hasDecay || hasVeto){
    PSPMTanode = event.GetSummary("pspmt:anode")->GetList();
    OtherAnode = event.GetSummary("generic:anode")->GetList();
  }

  for (auto it = OtherAnode.begin(); it!= OtherAnode.end(); it++){
    if (hasImplant){
      if ((*it)->GetChanID().HasTag("xa")){
        hi_xa = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("xb")){
	hi_xb = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("ya")){
	hi_ya = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("yb")){
	hi_yb = (*it)->GetCalibratedEnergy();
      }
    }
    if (hasDecay || hasVeto){
      if ((*it)->GetChanID().HasTag("xa")){
        low_xa = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("xb")){
	low_xb = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("ya")){
	low_ya = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("yb")){
	low_yb = (*it)->GetCalibratedEnergy();
      }

    }
  }

  for (auto it = PSPMTanode.begin(); it != PSPMTanode.end(); it++){
    if (hasImplant){
      if ((*it)->GetChanID().HasTag("xa")){
        low_xa = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("xb")){
	low_xb = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("ya")){
	low_ya = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("yb")){
	low_yb = (*it)->GetCalibratedEnergy();
      }
    }
    if (hasDecay || hasVeto){
      if ((*it)->GetChanID().HasTag("xa")){
        hi_xa = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("xb")){
	hi_xb = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("ya")){
	hi_ya = (*it)->GetCalibratedEnergy();
      }
      else if ((*it)->GetChanID().HasTag("yb")){
	hi_yb = (*it)->GetCalibratedEnergy();
      }

    }
    
    if ((*it)->GetChanID().HasTag("xa") && (*it)->GetCalibratedEnergy() > PspmtThreshold && /*(*it)->GetCalibratedEnergy() < 4000 &&*/ xa == 0){
      xa = (*it)->GetCalibratedEnergy();
      
    }
    if ((*it)->GetChanID().HasTag("xb") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&/* (*it)->GetCalibratedEnergy() < 4200 &&*/ xb == 0){
      xb = (*it)->GetCalibratedEnergy();

    }
    if ((*it)->GetChanID().HasTag("ya") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&/* (*it)->GetCalibratedEnergy() < 3800 &&*/ ya == 0){
      ya = (*it)->GetCalibratedEnergy();

    }
    if ((*it)->GetChanID().HasTag("yb") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&/* (*it)->GetCalibratedEnergy() < 4200 && */yb == 0){
      yb = (*it)->GetCalibratedEnergy();

    }
  }

  if (xa != 0 && xb != 0 && ya != 0 && yb != 0  && xa < 65000){
    position.first = (xa - xb) /(xa + xb);
    position.second = (ya - yb) / (ya + yb);
    if(hasImplant){
      pixel.first = Px * (5 * position.first + 1);
      pixel.second = Py * (5 * position.second + 1);
    }
    else if (hasDecay || hasVeto){
      pixel.first = Px * (5 * position.first + 1);
      pixel.second = Py * (5 * position.second + 1);
    }
    pos_x = position.first;
    pos_y = position.second;
    hasPosition = true;
  }
  else {
    if(hasDecay)
      plot(D_BAD_ANODE_MULT, PSPMTanode.size() + 10);
    if (hasImplant)
      plot(D_BAD_ANODE_MULT, PSPMTanode.size());
  }

  
  
  // PSPMT diagnostics


  // for implants
  if (hasPID && hasPosition){
    plot(D_DYNODE_MULT, dynodeLow.size());
    plot(D_ANODE_MULT, PSPMTanode.size());
    /*if ((*dynodeLow.begin())->GetCalibratedEnergy() > 5500 && 
      (*dynodeLow.begin())->GetCalibratedEnergy() < 5510){*/
      plot(DD_IMPLANT_POS, pspmtScale * position.first + pspmtOffset,
	   pspmtScale * position.second + pspmtOffset);
      // }
    if (pixel.first == pixel.second)
      plot(DD_IMPLANT_PIXELS, position.first * pspmtScale + pspmtOffset, position.second * pspmtScale + pspmtOffset);
    
    for(auto it = dynodeLow.begin(); it != dynodeLow.end(); it++){
      plot(D_IMPLANT_DYNODE, (*it)->GetCalibratedEnergy());
      plot(DD_ANODE_DYNODE, (xa+xb+ya+yb) / 4, (*it)->GetCalibratedEnergy());
      if (has72Co){
	plot(DD_XA_DYNODE, xa, (*it)->GetCalibratedEnergy());
	plot(DD_XB_DYNODE, xb, (*it)->GetCalibratedEnergy());
	plot(DD_YA_DYNODE, ya, (*it)->GetCalibratedEnergy());
	plot(DD_YB_DYNODE, yb, (*it)->GetCalibratedEnergy());
      }
    }
    
  }

  //for decays
  else if (hasDecay && hasPosition){
    plot(D_DYNODE_MULT, dynodeHi.size() + 10);
    plot(D_ANODE_MULT, PSPMTanode.size() + 10);
    plot(DD_DECAY_POS, pspmtScale * position.first + pspmtOffset, pspmtScale * position.second + pspmtOffset);
    if (pixel.first == pixel.second)
      plot(DD_DECAY_PIXELS, position.first * pspmtScale + pspmtOffset, position.second * pspmtScale + pspmtOffset);

    for(auto it = dynodeHi.begin(); it != dynodeHi.end(); it++){
      plot(D_DECAY_DYNODE, (*it)->GetCalibratedEnergy()/2);
      plot(DD_ANODE_DYNODE, (xa+xb+ya+yb) / 4, (*it)->GetCalibratedEnergy());
    }

  }

  //for events in the veto
  else if (hasVeto){
    plot(DD_VETO_POS, pspmtScale * position.first + pspmtOffset, pspmtScale * position.second + pspmtOffset);
    
    for (auto it = dynodeHi.begin(); it != dynodeHi.end(); it++){
      plot(D_VETO_DYNODE, (*it)->GetCalibratedEnergy());
    }
    
  }



  //---------------PSPMT position correlation---------------------

    static double pixel_time[2*Px][2*Py] = {};
    double decay_window = 79 * pow(10, 6) / 4;
    double decay_time = -1;

      int m = pixel.first -1; int n = pixel.second -1;

      if (m < 2 * Px && n < 2 * Py && m >= 0 && n >= 0){

        if (hasPID && hasPosition) {
            pixel_time[m][n] = (*dynodeLow.begin())->GetTimeSansCfd();
            //if (pixel.first == 12 && pixel.second == 12)


        }
        else if (hasDecay && hasPosition && m == 10 && n == 11){
            for(int x_it = 0; x_it < 2* Px; x_it++){
                for(int y_it = 0; y_it < 2 * Py; y_it++){
                    decay_time = (*dynodeHi.begin())->GetTimeSansCfd() - pixel_time[x_it][y_it];
                  if (decay_time > 0 && decay_time < decay_window)
                    plot(DD_HOTSPOT_MAP, x_it, y_it);
                }
            }
        }
      }


  


  //--------------- Ge plots -------------------------------------

  if (hasPID){
    for (auto it = geEvts.begin(); it != geEvts.end(); it++){
      plot(DD_PIN1_GE, (*it)->GetCalibratedEnergy(), pin1);
      plot(DD_GE_COR_TOF, i2pos1_cor_tof_pin1_i2n,
	   (*it)->GetCalibratedEnergy());
      plot(D_GE_IMPLANT, (*it)->GetCalibratedEnergy());
      if (has78Zn)
	plot(D_GE_ISOMER, (*it)->GetCalibratedEnergy()); 
    }
  }
  else if (hasDecay){
    for (auto it = geEvts.begin(); it != geEvts.end(); it++){
      plot(D_GE_DECAY, (*it)->GetCalibratedEnergy());
      if (has72CoDecay)
	plot(D_GE_CORR_DECAY, (*it)->GetCalibratedEnergy());
      //      else if(has71CoDecay)
      //	plot(D_GE_ANTI_CORR_DECAY, (*it)->GetCalibratedEnergy());
    }
  }

  //cout<<"3"<<endl;

  roottree->Fill();
  low_xa = 0;
  low_xb = 0;
  low_ya = 0;
  low_yb = 0;
  hi_xa = 0;
  hi_xb = 0;
  hi_ya = 0;
  hi_yb = 0;
  low_dynode = 0;
  hi_dynode = 0;
  timestamp = 0;
  
  //cout<<"4"<<endl;

  return(true);
}

double E14060Processor::CorrectTofByPos(const double &i2ns, const double &tof){

  double slope = -2.298;
  double intercept = 4056;

  return tof - slope * i2ns - intercept;
  
}
