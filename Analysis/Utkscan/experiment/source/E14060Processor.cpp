/** \file E14060Processor.cpp
 * \brief Experiment specific processor to handle the e14060 experiment at
 * the NSCL.
 *\author S. V. Paulauskas
 *\date September 15, 2016
 *\modified December 6, 2016 by A. Keeler
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

namespace dammIds {
    namespace experiment {
        const int DD_TRACE_MAX_DYNODE = 0;//!< Max value of the Dynode trace.
        const int DD_EPIN1_VS_TOF_PIN1_I2N = 1; //!<PIN1 vs. ToF between PIN1
//!< and I2N
        const int DD_EPIN1_VS_TOF_PIN1_I2S = 2; //!<PIN1 vs. ToF between PIN1
//!< and I2S
        const int DD_EPIN2_VS_TOF_PIN2_I2N = 3; //!<PIN2 vs. ToF between PIN2
//!< and I2N
        const int DD_EPIN2_VS_TOF_PIN2_I2S = 4; //!<PIN2 vs. ToF between PIN2
//!< and ISS
        const int DD_EPIN1_VS_TOF_I2N_I2S = 5; //!<PIN1 vs. ToF between I2N
//!< and I2S
  
        const int DD_TOF_I2NS_VS_TOF_PIN1_I2N = 6; //! ToF(I2N - I2S) vs. ToF(PIN1 - I2N)
        const int DD_PIN1_VS_PIN2 = 7; //!< PIN1 vs. PIN2

        const int DD_TOF_I2NS_VS_TOF_PIN1_I2N_CORR = 17; //! ToF(I2N - I2S) vs. ToF(PIN1 - I2N) - Corr
        const int DD_TOF_I2NS_PIN1_VS_TOF_PIN1_I2N = 25;
        const int DD_I2NS_VS_I2NS_PIN1 = 26;
        const int DD_I2NS_PIN1_VS_TOF_PIN1_I2N_CORR = 27;
        
        const int DD_TOF_I2NS_CORR_VS_E_PIN1 = 18;
        const int DD_TOF_I2NS_CORR_VS_E_PIN1_NOGATE = 24;

        const int DD_PIN1_VS_VETO = 21;
      
        namespace DECAY_GATED {
            const int DD_QDCVSTOF = 8; //!<QDC vs ToF
            const int DD_HAGRID = 9; //!< Location vs. HAGRiD En
            const int DD_NAI = 10; //!< Location vs. NaI En
	    const int DD_PSPMT_POS = 19; //!< PSPMT Position Gated with PIN
	    const int D_GE = 22; //Decay gated Ge
        }

        namespace PSPMT_GATED {
            const int DD_EPIN1_VS_TOF_PIN1_I2N = 11; //!<PIN1 vs. ToF between
//!< PIN1 and I2N anti-gated with the PSPMT signal
	    const int DD_EPIN1_VS_TOF_PIN1_I2N_CORR = 20; //!<PIN1 vs. ToF between
        }

        namespace IMPLANT_GATED {
            const int DD_QDCVSTOF = 12; //!<QDC vs ToF
            const int DD_EPIN1_VS_TOF_PIN1_I2N = 13; //!<PIN1 vs. ToF between
//!< PIN1 and I2N
        }

        namespace PIN_GATED {
            const int DD_PSPMT_POS = 14; //!< PSPMT Position Gated with PIN
        }

        namespace GE_GATED {
            const int DD_EPIN1_VS_TOF_PIN1_I2N = 15; //!<PIN1 vs. ToF between
//!< PIN1 and I2N
            const int DD_EPIN1_VS_GE = 16; //!< PIN1 Energy vs. Ge energy
	    const int DD_EPIN1_VS_TOF_PIN1_I2N_CORR = 23; //! PID Gated w/ GE
        }
    }
}//namespace dammIds



void E14060Processor::DeclarePlots(void) {
    DeclareHistogram2D(DD_TRACE_MAX_DYNODE, SC, S1, "Trace Max Dynode");
    DeclareHistogram2D(DD_EPIN1_VS_TOF_PIN1_I2N, SB, SB, "EPIN1 vs. ToF"
            "(I2N-PIN1)");
    DeclareHistogram2D(DD_EPIN1_VS_TOF_PIN1_I2S, SB, SB, "EPIN1 vs. ToF"
             "(I2S-PIN1)");
    DeclareHistogram2D(DD_EPIN2_VS_TOF_PIN2_I2N, SB, SB, "EPIN2 vs. ToF"
             "(I2N-PIN2)");
     DeclareHistogram2D(DD_EPIN2_VS_TOF_PIN2_I2S, SB, SB, "EPIN2 vs. ToF"
	     "(I2S-PIN2)");
    DeclareHistogram2D(DD_EPIN1_VS_TOF_I2N_I2S, SB, SB, "Si Energy vs. TOF"
	       "(I2N-I2S)");

    DeclareHistogram2D(DD_TOF_I2NS_VS_TOF_PIN1_I2N, SB, SB, "I2NS vs. TOF");
    DeclareHistogram2D(DD_TOF_I2NS_PIN1_VS_TOF_PIN1_I2N, SB, SB, "I2N-pin1 - I2S-pin1 vs. TOF");
    DeclareHistogram2D(DD_I2NS_VS_I2NS_PIN1, SB, SB, "I2NS vs. I2S-pin1 - I2N-pin1"); 
    DeclareHistogram2D(DD_TOF_I2NS_VS_TOF_PIN1_I2N_CORR, SB, SB, "I2NS vs. TOF CORR");
    DeclareHistogram2D(DD_I2NS_PIN1_VS_TOF_PIN1_I2N_CORR, SB, SB, "I2S-I2N pin1 vs. TOF CORR");
    DeclareHistogram2D(DD_TOF_I2NS_CORR_VS_E_PIN1, SB, SB, "DE1 vs. TOF CORR GATED");
    DeclareHistogram2D(DD_TOF_I2NS_CORR_VS_E_PIN1_NOGATE, SB, SB, "DE1 vs. TOF CORR NOT GATED");

    DeclareHistogram2D(DD_PIN1_VS_PIN2, SB, SB, "EPin1 vs. EPin2");

    //----- Histograms gated with decays
    DeclareHistogram2D(DECAY_GATED::DD_QDCVSTOF, SC, SD, "Decay - QDC vs. ToF");
    DeclareHistogram2D(DECAY_GATED::DD_PSPMT_POS, SB, SB,
                       "PSPMT Pos - Pin Gated");
    DeclareHistogram1D(DECAY_GATED::D_GE, SD, "Decay - Ge");
    DeclareHistogram2D(DECAY_GATED::DD_HAGRID, SC, S5, "Decay - HAGRiD");
    DeclareHistogram2D(DECAY_GATED::DD_NAI, SC, S5, "Decay - NaI");

    //----- Histograms gated on implants
    DeclareHistogram2D(IMPLANT_GATED::DD_QDCVSTOF, SC, SD, "Implant - QDC vs. "
            "ToF");
    DeclareHistogram2D(IMPLANT_GATED::DD_EPIN1_VS_TOF_PIN1_I2N, SB, SB,
                       "EPIN1 vs. ToF (I2N-PIN1) - YAP");

    //---------- Histograms Gated with the PIN
    DeclareHistogram2D(PIN_GATED::DD_PSPMT_POS, SB, SB,
                       "PSPMT Pos - Pin Gated");

    //----------- Histograms gated with the Clover
    DeclareHistogram2D(GE_GATED::DD_EPIN1_VS_TOF_PIN1_I2N, SB, SB,
                       "EPIN1 vs. ToF (I2N-PIN1) - YAP");
    DeclareHistogram2D(GE_GATED::DD_EPIN1_VS_GE, SC, SB, "EPIN1 vs. Ge");
    DeclareHistogram2D(GE_GATED::DD_EPIN1_VS_TOF_PIN1_I2N_CORR, SB, SB, "I2NS vs. TOF CORR");

    //----------- Histograms gated with the PSPMT
    DeclareHistogram2D(PSPMT_GATED::DD_EPIN1_VS_TOF_PIN1_I2N, SB, SB,
                       "EPIN1 vs. ToF (I2N-PIN1)- NO YAP");
    DeclareHistogram2D(PSPMT_GATED::DD_EPIN1_VS_TOF_PIN1_I2N_CORR, SB, SB,
                       "EPIN1 vs. ToF (I2N-PIN1) - CORR W YAP");
}

E14060Processor::E14060Processor(std::pair<double, double> &energyRange) :
        EventProcessor(OFFSET, RANGE, "E14060Processor") {
    SetAssociatedTypes();
    energyRange_ = energyRange;
}

void E14060Processor::SetAssociatedTypes() {
    associatedTypes.insert("vandle");
    associatedTypes.insert("hagrid");
    associatedTypes.insert("pspmt");
    associatedTypes.insert("ge");
    associatedTypes.insert("tac");
    associatedTypes.insert("si");
}

bool E14060Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return (false);

    static const double plotMult = 2;
    static const double plotOffset = 1000;
    static const double pspmtScale = 512; 
    static const double pspmtOffset = 512;

    map<string,double> pins_and_tacs;

    BarMap vbars;
    vector<ChanEvent *> geEvts;
    vector<vector<AddBackEvent>> geAddback;
    pair<double, double> position;
    pair<unsigned int, unsigned int> pixel;

    if (event.GetSummary("vandle")->GetList().size() != 0)
        vbars = ((VandleProcessor *) DetectorDriver::get()->
                GetProcessor("VandleProcessor"))->GetBars();
    if (event.GetSummary("pspmt:anode")->GetList().size() != 0) {
        position = ((PspmtProcessor *) DetectorDriver::get()->
                GetProcessor("PspmtProcessor"))->GetPosition("pixie");
        pixel = ((PspmtProcessor *) DetectorDriver::get()->
                GetProcessor("PspmtProcessor"))->GetPosition("pixie");
    }
    if (event.GetSummary("ge")->GetList().size() != 0) {
        geEvts = ((GeProcessor *) DetectorDriver::get()->
                GetProcessor("GeProcessor"))->GetGeEvents();
        geAddback = ((GeProcessor *) DetectorDriver::get()->
                GetProcessor("GeProcessor"))->GetAddbackEvents();
    }

    //--------- Obtain dynode information --------------------------------------
    static const vector<ChanEvent *> &dynode =
            event.GetSummary("pspmt:dynode")->GetList();
    static const vector<ChanEvent *> &dynodeClone = 
      event.GetSummary("generic:dynode")->GetList();

    TimingMapBuilder startbuilder(dynode);
    TimingMap tdynode = startbuilder.GetMap();

    //Loop over the dynode events to plot the maximum value
    for (TimingMap::const_iterator iterator3 = tdynode.begin();
         iterator3 != tdynode.end(); iterator3++)
      plot(DD_TRACE_MAX_DYNODE, iterator3->second.GetMaximumValue(), 0);

    for (vector<ChanEvent *>::const_iterator iterator1 = dynodeClone.begin();
	 iterator1 != dynodeClone.end(); iterator1++) {
      plot(DD_TRACE_MAX_DYNODE, (*iterator1)->GetTrace().GetMaxInfo().second, 1);
    }
    
    static const vector<ChanEvent *> &tac =
            event.GetSummary("tac", true)->GetList();
    static const vector<ChanEvent *> &pin =
	event.GetSummary("pin", true)->GetList();

    //Loop over the pin events in order to fill in the pins_and_tacs map
    for (vector<ChanEvent *>::const_iterator it = tac.begin(); it != tac.end(); it++) {
        pins_and_tacs.insert(make_pair((*it)->GetChanID().GetSubtype(),
                                       (*it)->GetCalibratedEnergy()));
        //cout<<(*it)->GetChanID().GetSubtype()<<endl<< (*it)->GetCalibratedEnergy()<<endl;
    }
    cout<<"Count:"<<pins_and_tacs.count("pin1_i2n")<<endl;
    //Loop over the pin events in order to fill in the pins_and_tacs map
    for (vector<ChanEvent *>::const_iterator it = pin.begin(); it != pin.end(); it++)
	pins_and_tacs.insert(make_pair((*it)->GetChanID().GetSubtype(), 
				       (*it)->GetCalibratedEnergy()));

    //Basic correlation information
    bool hasIon = pin.size() != 0;
    bool hasVeto = event.GetSummary("generic:veto")->GetMult() != 0;
    bool hasImplant = hasIon && dynode.size() != 0;
    bool hasDecay = !hasIon && dynode.size() != 0;
    bool hasStrictDecay = hasDecay && !hasVeto;

    //--------Obtain information about the pins and tacs of interest
    //double pin1_i2n = FindPinOrTacEnergy(pins_and_tacs, "pin1_i2n");
    double pin1_i2n = pins_and_tacs.find("pin1_i2n")->second;
    //cout<<"Found:"<<pin1_i2n<<endl;
    //We put this here to gate out the overflow bins for the TAC
    if(pin1_i2n > 2000 && pin1_i2n < 100)
	pin1_i2n = 0.0;
    double pin1_i2s = FindPinOrTacEnergy(pins_and_tacs, "pin1_i2s");
    if(pin1_i2s > 2000 && pin1_i2s < 100)
	pin1_i2s = 0.0;
    double pin1 = FindPinOrTacEnergy(pins_and_tacs, "de1");
    double pin2 = FindPinOrTacEnergy(pins_and_tacs, "de2");
    double i2ns = FindPinOrTacEnergy(pins_and_tacs, "i2n_i2s");
    double i2ns_pin1 = 1.8 * (pin1_i2s - pin1_i2n);
    bool isCentralI2NS = i2ns > 1200 && i2ns < 1450; //Gate on the i2ns tac
    double i2ns_cor_tof = 0.0;
    if(pin1_i2n != 0.0 && i2ns != 0.0)
	i2ns_cor_tof = CorrectToFByI2Pos("pin1_i2n", pin1_i2n, i2ns);
    double i2ns_pin1_cor_tof = 0.0;
    if(pin1_i2n != 0.0 && i2ns_pin1 != 0.0)
	i2ns_pin1_cor_tof = CorrectToFByI2Pos("pin1_i2n", pin1_i2n, i2ns_pin1);
    bool is72co = histo.BananaTest(8, i2ns_cor_tof, pin1);
    bool is78zn = histo.BananaTest(3, i2ns_cor_tof, pin1);

    //-------------------- GE INFORMATION --------------------------------------
    //Loop over the Ge events to ensure that we had ourselves something in
    // the region of interest. We also take the time to plot the PIN1 energy
    // vs. Ge Energy
    bool hasGe = false;
    for (vector<ChanEvent *>::const_iterator iterator2 = geEvts.begin();
         iterator2 != geEvts.end(); iterator2++) {
        if ((*iterator2)->GetCalibratedEnergy() > energyRange_.first &&
            (*iterator2)->GetCalibratedEnergy() < energyRange_.second) {
            hasGe = true;
	}
	
	//Plots the decay gated gamma ray spectrum
	if(is78zn)
	    plot(DECAY_GATED::D_GE, (*iterator2)->GetCalibratedEnergy());
	plot(GE_GATED::DD_EPIN1_VS_GE, (*iterator2)->GetCalibratedEnergy(), pin1);
    }
    //-------------------- PLOTTING PIN ENERGIES -------------------------------
    double delta = 0.87558748 * pin2 + 8.72618557 - pin1;

   if (pin1 != 0 && pin2 != 0 && delta < 40 && delta > -40)
        plot(DD_PIN1_VS_PIN2, pin2, pin1);

 


    //------------------------ PLOTTING PID ------------------------------------
    if (delta < 40 && delta > -40){
      if (pin1 != 0) {
	if(pin1_i2n != 0 ) {
        cout<<pin1_i2n<<endl;
	  plot(DD_EPIN1_VS_TOF_PIN1_I2N, pin1_i2n, pin1);
	  if (hasGe)
	    plot(GE_GATED::DD_EPIN1_VS_TOF_PIN1_I2N, pin1_i2n, pin1);
	  if (hasImplant)
	    plot(IMPLANT_GATED::DD_EPIN1_VS_TOF_PIN1_I2N, pin1_i2n, pin1);
	  if(!hasImplant)
	    plot(PSPMT_GATED::DD_EPIN1_VS_TOF_PIN1_I2N, pin1_i2n, pin1);
	  
	  //Spectrum to use to give the position correction for the I2N/S TAC. 
	  plot(DD_TOF_I2NS_VS_TOF_PIN1_I2N, pin1_i2n, i2ns);
	  plot(DD_TOF_I2NS_PIN1_VS_TOF_PIN1_I2N, pin1_i2n,i2ns_pin1);
	  plot(DD_I2NS_VS_I2NS_PIN1, i2ns_pin1, i2ns);
	  plot(DD_I2NS_PIN1_VS_TOF_PIN1_I2N_CORR, i2ns_pin1_cor_tof, i2ns_pin1);

	  //PID for the corrected I2NS Position without any gates on the position
	  plot(DD_TOF_I2NS_CORR_VS_E_PIN1_NOGATE, i2ns_cor_tof, pin1);
	    //Gating on the TAC for I2N/S to remove wings and over/under flows
	  if( isCentralI2NS && i2ns_cor_tof != 0.0 ) {
	    plot(DD_TOF_I2NS_VS_TOF_PIN1_I2N_CORR, i2ns_cor_tof, i2ns);
	    plot(DD_TOF_I2NS_CORR_VS_E_PIN1, i2ns_cor_tof, pin1);
	    if (hasImplant) 
	      plot(PSPMT_GATED::DD_EPIN1_VS_TOF_PIN1_I2N_CORR, i2ns_cor_tof, pin1);
	    if(hasGe)
	      plot(GE_GATED::DD_EPIN1_VS_TOF_PIN1_I2N_CORR, i2ns_cor_tof, pin1);
	  }
	} //pin1_i2n != 0
	if (FindPinOrTacEnergy(pins_and_tacs,"pin1_i2s") != 0.0)
	  plot(DD_EPIN1_VS_TOF_PIN1_I2S, FindPinOrTacEnergy(pins_and_tacs,"pin1_i2s"), pin1);
      }//pin1 != 0
      
      //Spectra for Pin2
      if(FindPinOrTacEnergy(pins_and_tacs,"de2") != 0.0) {
	if (FindPinOrTacEnergy(pins_and_tacs,"pin2_i2n") != 0.0)
	  plot(DD_EPIN2_VS_TOF_PIN2_I2N, FindPinOrTacEnergy(pins_and_tacs,"pin2_i2n"), 
	       FindPinOrTacEnergy(pins_and_tacs,"de2"));
	if (FindPinOrTacEnergy(pins_and_tacs,"pin2_i2s") != 0.0)
	  plot(DD_EPIN2_VS_TOF_PIN2_I2S, FindPinOrTacEnergy(pins_and_tacs,"pin2_i2s"), 
	       FindPinOrTacEnergy(pins_and_tacs,"de2"));
      }
    }
    // ----------------- PLOTTING POSITION OF PSPMT ----------------------------
    if(hasImplant && position.first != 0.0 && position.second != 0.0)
      plot(PIN_GATED::DD_PSPMT_POS, position.first * pspmtScale + pspmtOffset, 
	   position.second * pspmtScale + pspmtOffset);
    if(hasStrictDecay && position.first != 0.0 && position.second != 0.0)
      plot(DECAY_GATED::DD_PSPMT_POS, position.first * pspmtScale + pspmtOffset, 
	   position.second * pspmtScale + pspmtOffset);

    // ----------------- HAGRID AND NAI ----------------------------------------
    static const vector<ChanEvent *> &hagrid =
            event.GetSummary("hagrid")->GetList();

    static const vector<ChanEvent *> &nai =
            event.GetSummary("nai")->GetList();

    if (hasDecay) {
        for (vector<ChanEvent *>::const_iterator iterator1 = hagrid.begin();
             iterator1 != hagrid.end(); iterator1++)
            plot(DECAY_GATED::DD_HAGRID, (*iterator1)->GetCalibratedEnergy(),
                 (*iterator1)->GetChanID().GetLocation());
        for (vector<ChanEvent *>::const_iterator iterator1 = nai.begin();
             iterator1 != nai.end(); iterator1++)
            plot(DECAY_GATED::DD_NAI, (*iterator1)->GetCalibratedEnergy(),
                 (*iterator1)->GetChanID().GetLocation());
    }

    //--------------- VANDLE ---------------------------------------------------
    for (BarMap::iterator it = vbars.begin(); it != vbars.end(); it++) {
        BarDetector bar = (*it).second;

        if (!bar.GetHasEvent())
            continue;

        TimingCalibration cal = bar.GetCalibration();

        for (TimingMap::iterator itStart = tdynode.begin();
             itStart != tdynode.end(); itStart++) {

            if (!(*itStart).second.GetIsValid())
                continue;

            unsigned int startLoc = (*itStart).first.first;
            HighResTimingData start = (*itStart).second;

            double tof = bar.GetCorTimeAve() -
                         start.GetWalkCorrectedTime() + cal.GetTofOffset(startLoc);

            if (hasDecay)
                plot(DECAY_GATED::DD_QDCVSTOF, tof * plotMult + plotOffset,
                     bar.GetQdc());
            if (hasImplant)
                plot(IMPLANT_GATED::DD_QDCVSTOF, tof * plotMult + plotOffset,
                     bar.GetQdc());
        }//for(TimingMap::iterator start = tdynode.begin();
    } //for(BarMap::iterator it = vbars.begin()
    
    
    EndProcess();
    return (true);
}

double E14060Processor::FindPinOrTacEnergy(const std::map<string,double> &mp,
					   const std::string &key) {
    map<string,double>::const_iterator it = mp.find(key);
    //cout<<"key="<<key<<endl<<"Energy:"<<(*it).second<<endl;
    if(it == mp.end())
	return 0.0;
    else
	return (*it).second;
}

double E14060Processor::CorrectToFByI2Pos(const std::string &name, 
					  const double &tof,
					  const double &i2ns) {
    double slope = 0.0;
    double intercept = 0.0;
    double dammOffset = 1000.;
    if(name == "pin1_i2n") {
      	slope = -2.298;
      //	intercept = 4684.7;
      //  slope =-0.448;
        intercept = 5056;
    } else if (name == "pin1_i2s") {
    	slope = 1.;
    	intercept = 0.0;
    } else if (name == "pin2_i2n") {
    	slope = 1.;
    	intercept = 0.0;
    } else if (name == "pin2_i2s") {
    	slope = 1.;
    	intercept = 0.0;
    }
    return tof - slope * i2ns - intercept + dammOffset;
}
