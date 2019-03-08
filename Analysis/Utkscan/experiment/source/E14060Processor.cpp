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
        const int D_DYNODE_DELAY = 19;
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

    DeclareHistogram1D(D_DYNODE_DELAY, S9, "Time difference for multiple dynode events");
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
    DeclareHistogram2D(DD_HOTSPOT_MAP, SA, SA, "Implants correlated with single-pixel decays");
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
    hi_dynode = 0;
    hi_dynode_time = 0;
    hi_dynode_mult = 0;
    low_dynode = 0;
    low_dynode_time = 0;
    low_dynode_mult = 0;
    low_dynode_tr_max = 0;
    low_dynode_tr_max = {};
    

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
    roottree->Branch("hi_dynode", &hi_dynode);
    roottree->Branch("hi_dynode_time", &hi_dynode_time);
    roottree->Branch("hi_dynode_mult", &hi_dynode_mult);
    roottree->Branch("low_dynode", &low_dynode);
    roottree->Branch("low_dynode_time", &low_dynode_time);
    roottree->Branch("low_dynode_mult", &low_dynode_mult);
    roottree->Branch("low_dynode_tr_max", &low_dynode_tr_max);
    roottree->Branch("low_dynode_trace", &low_dynode_trace);
    roottree->Branch("event_type", &eventType);
    roottree->Branch("PID", &pid_event);
    roottree->Branch("Events", &current_event);
    roottree->Branch("past_events", &pastEvents);


}

void E14060Processor::SetAssociatedTypes() {
  //associatedTypes.insert("vandle");
  //associatedTypes.insert("hagrid");
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
    event.GetSummary("pspmt:dynode_low")->GetList();
  static const vector<ChanEvent *> &dynodeHi =
    event.GetSummary("pspmt:dynode_high")->GetList();

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
    bool hasDynodeHi = !dynodeHi.empty() && (*dynodeHi.begin())->GetCalibratedEnergy() > 1 && (*dynodeHi.begin())->GetCalibratedEnergy() < 16000;

    bool hasImplantReject = hasLightIon || hasVeto;
    bool hasDecayReject = hasVeto || hasIon;

    bool hasImplant = hasIon && hasDynodeLow && !hasImplantReject;
    bool hasDecay = hasDynodeHi && !hasDecayReject;
    if(hasImplant)
        timestamp = (*dynodeLow.begin())->GetTimeSansCfd();
    if(hasDecay)
        timestamp = (*dynodeHi.begin())->GetTimeSansCfd();


    if (hasImplant)
        eventType = "implant";
    else if (hasDecay)
        eventType = "decay";
    else if (hasVeto)
        eventType = "veto";
    else
        eventType = " ";

    bool has72Co = false;
    bool has71Co = false;
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
                if(i2pos1_cor_tof_pin1_i2n < 773 && i2pos1_cor_tof_pin1_i2n > 623){
                   if(pin1 < 531 && pin1 > 493)
                       has71Co = true;
            }
            if (i2pos1_cor_tof_pin1_i2n < 1016 && i2pos1_cor_tof_pin1_i2n > 885){
                if (pin1 < 664 && pin1 > 617)
                   has78Zn = true;
            }
        }
    }

    if(hasPID){
      pid_event.pid = true;
      pid_event.pin = pin1;
      pid_event.beam_tof = pin1_i2n;
      pid_event.i2pos = i2pos1;
      pid_event.corrected_beam_tof = i2pos1_cor_tof_pin1_i2n;
    } else{
      pid_event = default_pid;
    }



  //-------------------Calculate positions for Pspmt---------------

  double xa = 0; double xb = 0; double ya = 0; double yb = 0;
  //const int Px = 12; const int Py = 12;
  double PspmtThreshold = 50;
  bool hasPosition = false;

    vector<ChanEvent *> PSPMTanodeHi;
    PSPMTanodeHi = event.GetSummary("pspmt:anode_high")->GetList();
    vector<ChanEvent *> PSPMTanodeLow;
    PSPMTanodeLow = event.GetSummary("pspmt:anode_low")->GetList();


  for (auto it = PSPMTanodeHi.begin(); it!= PSPMTanodeHi.end(); it++){
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
      if(hasDecay || hasVeto){
          if ((*it)->GetChanID().HasTag("xa") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&
                (*it)->GetCalibratedEnergy() < 65000 && xa == 0){
              xa = (*it)->GetCalibratedEnergy();

          }
          if ((*it)->GetChanID().HasTag("xb") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&
                (*it)->GetCalibratedEnergy() < 65000 && xb == 0){
              xb = (*it)->GetCalibratedEnergy();

          }
          if ((*it)->GetChanID().HasTag("ya") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&
                (*it)->GetCalibratedEnergy() < 65000 && ya == 0){
              ya = (*it)->GetCalibratedEnergy();

          }
          if ((*it)->GetChanID().HasTag("yb") && (*it)->GetCalibratedEnergy() > PspmtThreshold &&
                (*it)->GetCalibratedEnergy() < 65000 && yb == 0){
              yb = (*it)->GetCalibratedEnergy();

          }
      }

  }

  for (auto it = PSPMTanodeLow.begin(); it != PSPMTanodeLow.end(); it++){
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

      if(hasImplant){
        if ((*it)->GetChanID().HasTag("xa") && (*it)->GetCalibratedEnergy() > PspmtThreshold && xa == 0){
            xa = (*it)->GetCalibratedEnergy();

        }
        if ((*it)->GetChanID().HasTag("xb") && (*it)->GetCalibratedEnergy() > PspmtThreshold && xb == 0){
            xb = (*it)->GetCalibratedEnergy();

        }
        if ((*it)->GetChanID().HasTag("ya") && (*it)->GetCalibratedEnergy() > PspmtThreshold && ya == 0){
            ya = (*it)->GetCalibratedEnergy();

        }
        if ((*it)->GetChanID().HasTag("yb") && (*it)->GetCalibratedEnergy() > PspmtThreshold && yb == 0){
            yb = (*it)->GetCalibratedEnergy();

        }
    }
  }

  if (xa != 0 && xb != 0 && ya != 0 && yb != 0){

    if(hasImplant){
        position.first = ((xa - xb + 0.0005 * (pow(xa,2) - pow(xb,2)))
	    / (xa + xb + 0.0005 * (pow(xa,2) + pow(xb, 2)))) * 0.87 + 0.014;
        position.second = ((ya - yb + 0.0006 * (pow(ya,2) - pow(yb,2)))
	    / (ya + yb + 0.0006 * (pow(ya,2) + pow(yb, 2)))) * 0.85 - 0.007;
        pixel.first = Px * (3.33 * position.first + 1);
        pixel.second = Py * (3.33 * position.second + 1);
    }
    else if (hasDecay || hasVeto){
        position.first = (xa - xb) /(xa + xb);
        position.second = (ya - yb) / (ya + yb);
        pixel.first = Px * (3.33 * position.first + 1);
        pixel.second = Py * (3.33 * position.second + 1);
    }
    hasPosition = true;
  }


  
  
  // PSPMT diagnostics

  // for implants
    if (hasPID && hasPosition){
        plot(D_DYNODE_MULT, dynodeLow.size());
      //plot(D_ANODE_MULT, PSPMTanode.size());

        plot(DD_IMPLANT_POS, pspmtScale * position.first + pspmtOffset,
           pspmtScale * position.second + pspmtOffset);

      if (dynodeLow.size() > 1)
          plot(D_DYNODE_DELAY,(dynodeLow.back())->GetTimeSansCfd() -
                (*dynodeLow.begin())->GetTimeSansCfd());

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
    //plot(D_ANODE_MULT, PSPMTanode.size() + 10);
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
    //PspmtEvent current_event;



    if (hasPosition) {


        current_event.implant = hasPID;
        current_event.decay = hasDecay;

        current_event.pidEvent = pid_event;

        /*       pid_event.pid = true;
        pid_event.pin = pin1;
        pid_event.beam_tof = pin1_i2n;
        pid_event.i2pos = i2pos1;
        pid_event.corrected_beam_tof = i2pos1_cor_tof_pin1_i2n;
        */
        current_event.x_position = position.first;
        current_event.y_position = position.second;
        current_event.x_pixel = pixel.first;
        current_event.y_pixel = pixel.second;
        current_event.event_time = timestamp;
        current_event.pixel_num = current_event.x_pixel * 24 + current_event.y_pixel;
    } else {
        current_event = defaultStruct;
    }
        if(dynodeLow.size() > 0){
            low_dynode = (*dynodeLow.begin())->GetCalibratedEnergy();
            low_dynode_time = (*dynodeLow.begin())->GetTimeSansCfd();
            low_dynode_mult = dynodeLow.size();
            low_dynode_tr_max = dynodeLow.front()->GetTrace().GetMaxInfo().second;
            //low_dynode_trace = dynodeLow.front()->GetTrace();
        }
        if(dynodeHi.size() > 0){
            hi_dynode = (*dynodeHi.begin())->GetCalibratedEnergy();
            hi_dynode_time = (*dynodeHi.begin())->GetTimeSansCfd();
            hi_dynode_mult = dynodeHi.size();
        }




    //--------------- Ge plots -------------------------------------

    for (auto it = geEvts.begin(); it != geEvts.end(); it++){
      current_event.gammaEvents.emplace_back((*it)->GetTimeSansCfd(),(*it)->GetCalibratedEnergy());
      if(hasPID){
        plot(DD_PIN1_GE, (*it)->GetCalibratedEnergy(), pin1);
        plot(DD_GE_COR_TOF, i2pos1_cor_tof_pin1_i2n, (*it)->GetCalibratedEnergy());
        plot(D_GE_IMPLANT, (*it)->GetCalibratedEnergy());
        if (has78Zn)
          plot(D_GE_ISOMER, (*it)->GetCalibratedEnergy());
      }

      else if (hasDecay){
        plot(D_GE_DECAY, (*it)->GetCalibratedEnergy());
        //if (has72CoDecay)
        //plot(D_GE_CORR_DECAY, (*it)->GetCalibratedEnergy());
        //      else if(has71CoDecay)
        //	plot(D_GE_ANTI_CORR_DECAY, (*it)->GetCalibratedEnergy());
      }
    }



  //---------------PSPMT position correlation---------------------


    static std::vector<PspmtEvent> past_events;
    double decay_window = 79 * pow(10, 6) / 4;

    if(hasPosition && current_event.implant)
        past_events.emplace_back(current_event);

    if (hasPosition && current_event.decay){
        for(auto it = past_events.begin(); it != past_events.end(); it++) {
            if (past_events.empty()) {
                PspmtEvent tmp;
                past_events.emplace_back(tmp);
                break;
            }

            if (current_event.event_time - (*it).event_time < decay_window) {
                plot(DD_HOTSPOT_MAP, ((*it).x_position * pspmtScale + pspmtOffset),
                     (*it).y_position * pspmtScale + pspmtOffset);
                    pastEvents.emplace_back(*it);

            }

            if (current_event.event_time - (*it).event_time > decay_window)
                past_events.erase(it);

        }
    }



    roottree->Fill();
    low_xa = 0;
    low_xb = 0;
    low_ya = 0;
    low_yb = 0;
    hi_xa = 0;
    hi_xb = 0;
    hi_ya = 0;
    hi_yb = 0;
    hi_dynode = 0;
    hi_dynode_time = 0;
    hi_dynode_mult = 0;
    low_dynode = 0;
    low_dynode_time = 0;
    low_dynode_mult = 0;
    low_dynode_tr_max = 0;
    low_dynode_trace.clear();

    pastEvents.clear();
    //gammaEvents.clear()
  

  return(true);
}

double E14060Processor::CorrectTofByPos(const double &i2ns, const double &tof){

  double slope = -2.298;
  double intercept = 4056;

  return tof - slope * i2ns - intercept;
  
}


