/** \file E14060Processor.cpp
 * \cleaned-up Experiment-specific processor for the e14060 experiment at 
 * the NSCL.  Based on a processor written by S. V. Paulauskas.
 *\ author A. Keeler
 *\ date April 3, 2017
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
        const int DD_MAX_DYNODE_TRACE = 0;  //!<Max value of the Dynode trace

        // Pin1 - Pin2 comparisons
        const int DD_PIN1_VS_PIN2 = 1;
        const int DD_PIN1_VS_PIN2_CUT = 2;

        // PID and position correction plots
        const int DD_EPIN1_VS_TOF_PIN1_I2N = 3;
        const int DD_TOF_PIN1_I2N_VS_I2NS = 4;
        const int DD_TOF_PIN1_I2N_VS_I2S_I2N = 5;
        const int DD_TOF_PIN1_I2N_CORR_VS_I2NS = 6;
        const int DD_TOF_PIN1_I2N_CORR_VS_I2S_I2N = 7;
        const int DD_EPIN1_VS_TOF_PIN1_I2NS_CORR = 8;
        const int DD_EPIN1_VS_TOF_PIN1_I2S_I2N_CORR = 9;

        // PSPMT position spectra
        const int DD_IMPLANT_PSPMT = 10;
        const int DD_ION_PSPMT = 11;
        const int DD_DECAY_PSPMT = 12;

        // Ge Information
        const int DD_EPIN1_VS_GE = 13;
        const int D_GE_DECAY_GATED = 14;
        const int D_GE_CORR_DECAY = 15;
        const int D_GE_SELECT_GATED = 16;

        // Veto Information
        const int D_VETO = 18;
        const int D_HAS_VETO = 19;

        // Correlation plots
        const int DD_CORR_DECAY_PSPMT = 20;
        const int DD_CORR_DECAY_PIXEL = 21;
        const int D_DECAY_DIST = 22;
        const int D_DECAY_TIME = 23;
    }
}//namespace dammIds

void E14060Processor::DeclarePlots(void) {
    DeclareHistogram2D(DD_MAX_DYNODE_TRACE, SC, S1, "Trace MAX Dynode");
    DeclareHistogram2D(DD_PIN1_VS_PIN2, SB, SB, "dE1 vs. dE2");
    DeclareHistogram2D(DD_PIN1_VS_PIN2_CUT, SB, SB, "dE1 vs. dE2 - cut");
    DeclareHistogram2D(DD_EPIN1_VS_TOF_PIN1_I2N, SB, SB, "dE1 vs TOF PIN1-I2n");
    DeclareHistogram2D(DD_TOF_PIN1_I2N_VS_I2NS, SB, SB, "TOF PIN1-I2n vs. I2ns position");
    DeclareHistogram2D(DD_TOF_PIN1_I2N_VS_I2S_I2N, SB, SB, "TOF PIN1-I2n vs. I2s-I2n position");
    DeclareHistogram2D(DD_TOF_PIN1_I2N_CORR_VS_I2NS, SB, SB, "corrected TOF PIN1-I2n vs. I2ns position");
    DeclareHistogram2D(DD_TOF_PIN1_I2N_CORR_VS_I2S_I2N, SB, SB, "corrected TOF PIN1-I2n vs. I2s-I2n position");
    DeclareHistogram2D(DD_EPIN1_VS_TOF_PIN1_I2NS_CORR, SB, SB, "dE1 vs. I2NS corrected TOF PIN1-I2n");
    DeclareHistogram2D(DD_EPIN1_VS_TOF_PIN1_I2S_I2N_CORR, SB, SB, "dE1 vs. I2S-I2N corrected TOF PIN1-I2n");

    DeclareHistogram2D(DD_IMPLANT_PSPMT, SB, SB, "Implant positions (PSPMT gated by PINs)");
    DeclareHistogram2D(DD_ION_PSPMT, SB, SB, "Light Ion positions (PSPMT gated on Veto)");
    DeclareHistogram2D(DD_DECAY_PSPMT, SB, SB, "Decay positions (PSPMT gated on decays)");

    DeclareHistogram2D(DD_EPIN1_VS_GE, SC, SB, "dE1 vs. Ge");
    DeclareHistogram1D(D_GE_DECAY_GATED, SD, "Ge - Ungated");
    DeclareHistogram1D(D_GE_CORR_DECAY, SD, "Ge - Correlated Decays");
    DeclareHistogram1D(D_GE_SELECT_GATED, SD, "Ge - 72Co decay-gated");

    DeclareHistogram1D(D_VETO, SE, "Veto events");
    DeclareHistogram1D(D_HAS_VETO, SE, "'true' Veto Events");

    DeclareHistogram2D(DD_CORR_DECAY_PSPMT, SA, SA, "Decay positions Correlated with Implants");
    DeclareHistogram2D(DD_CORR_DECAY_PIXEL, S5, S5, "Decay pixels Correlated with Implants");
    DeclareHistogram1D(D_DECAY_DIST, S8, "Distance between implant and decay");
    DeclareHistogram1D(D_DECAY_TIME, SD, "Time Between implant and decay");
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
    static const double pspmtScale = 500;
    static const double pspmtOffset = 500;

    map<string,double> pins_and_tacs;

    BarMap vbars;
    vector<ChanEvent *> geEvts;
    vector<ChanEvent *> vetoEvts;
    pair <double, double> position;
    pair <int, int> pixel;
    static map<pair<int, int>, pair<double, double>> implant_map;

    if (event.GetSummary("vandle")->GetList().size() != 0)
	vbars = ((VandleProcessor *) DetectorDriver::get()->
		 GetProcessor("VandleProcessor"))->GetBars();
    if (event.GetSummary("pspmt:anode")->GetList().size() != 0) {
        position = ((PspmtProcessor *) DetectorDriver::get()->
		 GetProcessor("PspmtProcessor"))->GetPosition("pixie");
        pixel = ((PspmtProcessor *) DetectorDriver::get()->
		 GetProcessor("PspmtProcessor"))->GetPixel("pixie");
    }

    if (event.GetSummary("ge")->GetList().size() != 0) {
        //static const vector<ChanEvent *> geEvts = ((GeProcessor *) DetectorDriver::get()->GetProcessor("GeProcessor"))->GetGeEvents();
        //static const vector<ChanEvent *> &
        geEvts =event.GetSummary("ge")->GetList();
    }

    if (event.GetSummary("generic:veto")->GetList().size() != 0) {
        //static const vector<ChanEvent *> &
        vetoEvts = event.GetSummary("generic:veto")->GetList();
    }

    //-------------- Obtain Dynode Information ----------------------------
    static const vector<ChanEvent *> &dynode = 
      event.GetSummary("pspmt:dynode")->GetList();
    static const vector<ChanEvent *> &dynodeClone =
      event.GetSummary("generic:dynode")->GetList();

    TimingMapBuilder startbuilder(dynode);
    TimingMap tdynode = startbuilder.GetMap();

    //Looping over the dynode events to plot the maximum value
    for (TimingMap::const_iterator iterator1 = tdynode.begin();
	 iterator1 != tdynode.end(); iterator1++)
      plot(DD_MAX_DYNODE_TRACE, iterator1->second.GetMaximumValue(), 0);
    for (vector<ChanEvent *>::const_iterator iterator2 = dynodeClone.begin();
	 iterator2 != dynodeClone.end(); iterator2++) {
      plot(DD_MAX_DYNODE_TRACE, (*iterator2)->GetTrace().GetMaxInfo().second, 1);
    }

    static const vector<ChanEvent *> &tac =
      event.GetSummary("tac", true)->GetList();
    static const vector<ChanEvent *> &pin = 
      event.GetSummary("pin", true)-> GetList();

    //Loop over the tac events in order to fill in the pins_and_tacs map
    for (vector<ChanEvent *>::const_iterator it = tac.begin(); it != tac.end(); 
	 it++) {
      pins_and_tacs.insert(make_pair((*it)->GetChanID().GetSubtype(),
				     (*it)->GetCalibratedEnergy()));
    }

    //Loop over the pin events in order to fill in the pins_and_tacs map
    for (vector<ChanEvent *>::const_iterator it = pin.begin(); it != pin.end();
	 it++)
      pins_and_tacs.insert(make_pair((*it)->GetChanID().GetSubtype(),
				    (*it)->GetCalibratedEnergy()));

    //------------ Obtain information about the pins or tacs of interest
    double pin1_i2n = pins_and_tacs.find("pin1_i2n")->second;
    double pin1_i2s = pins_and_tacs.find("pin1_i2s")->second;
    double pin2_i2n = pins_and_tacs.find("pin2_i2n")->second;
    double pin2_i2s = pins_and_tacs.find("pin2_i2s")->second;

    double i2ns = pins_and_tacs.find("i2n_i2s")->second;
    double i2pos1 = 1.8 * (pin1_i2s - pin1_i2n);
    double i2pos2 = 1.8 * (pin2_i2s - pin2_i2n);

    double pin1 = pins_and_tacs.find("de1")->second;
    double pin2 = pins_and_tacs.find("de2")->second;
    
    //All corrected ToFs are from I2N to one of the PINs
    double pin1_i2ns_cor_tof = 0.0;
    if (pin1_i2n != 0 && i2ns != 0)
      pin1_i2ns_cor_tof = CorrectToFByI2Pos("pin1_i2n", pin1_i2n, i2ns);
    double pin1_i2pos1_cor_tof = 0.0;
    if (pin1_i2n != 0 && i2pos1 != 0)
      pin1_i2pos1_cor_tof = CorrectToFByI2Pos("pin1_i2n", pin1_i2n, i2pos1);
    double pin2_i2pos2_cor_tof = 0.0;
    if (pin2_i2n != 0 && i2pos2 != 0)
      pin2_i2pos2_cor_tof = CorrectToFByI2Pos("pin2_i2n", pin2_i2n, i2pos2);

    // Basic Correlation information
    bool hasIon = pin.size() != 0;
    bool hasVeto = false;
    vector<ChanEvent *>::const_iterator it = vetoEvts.begin();
    while (!hasVeto && it != vetoEvts.end()){
        if ((*it)->GetCalibratedEnergy() > 400)
            hasVeto = true;
        it++;
    }
    //bool hasVeto = event.GetSummary("generic:veto")->GetMult() != 0;
    bool hasImplant = hasIon && dynodeClone.size() != 0 && !hasVeto;
    bool hasDecay = !hasIon && dynode.size() != 0 && !hasVeto;
    bool has72Co = false;

    //------------------ Plotting Ge Information -----------------------

    for (vector<ChanEvent *>::const_iterator it = vetoEvts.begin(); it != vetoEvts.end(); it++){
        plot(D_VETO, (*it)->GetCalibratedEnergy());
        if (hasVeto)
            plot(D_HAS_VETO, (*it)->GetCalibratedEnergy());
    }

    bool hasGe = false;

    for (vector<ChanEvent *>::const_iterator iterator2 = geEvts.begin();
         iterator2 != geEvts.end(); iterator2++) {
        if ((*iterator2)->GetCalibratedEnergy() > energyRange_.first &&
            (*iterator2)->GetCalibratedEnergy() < energyRange_.second) {
            hasGe = true;
        }
        plot(DD_EPIN1_VS_GE, (*iterator2)->GetCalibratedEnergy(), pin1);
    }


    //------------------ PLOTTING PIN ENERGIES -------------------------
    
    plot(DD_PIN1_VS_PIN2, pin2, pin1);
    double delta = 0.87558748 * pin2 + 8.72618557 - pin1;

    if (pin1 != 0 && pin2 != 0 && delta < 40 && delta > -40)
      plot(DD_PIN1_VS_PIN2_CUT, pin2, pin1);


    //---------------------- PLOTTING PID ------------------------------
    if (delta < 40 && delta > -40) {
        if (pin1 != 0 && pin1_i2n != 0) {
            //Plot PID and spectra for applying position correction
            plot(DD_EPIN1_VS_TOF_PIN1_I2N, pin1_i2n, pin1);
            //Plot i2ns position correction
            plot(DD_TOF_PIN1_I2N_VS_I2NS, i2ns, pin1_i2n);
            plot(DD_TOF_PIN1_I2N_CORR_VS_I2NS, i2ns, pin1_i2ns_cor_tof);
            //Plot Tof difference position correction
            plot(DD_TOF_PIN1_I2N_VS_I2S_I2N, i2pos1, pin1_i2n);
            plot(DD_TOF_PIN1_I2N_CORR_VS_I2S_I2N, i2pos1, pin1_i2pos1_cor_tof);
            //Plot corrected PIDs
            plot(DD_EPIN1_VS_TOF_PIN1_I2NS_CORR, pin1_i2ns_cor_tof, pin1);
            plot(DD_EPIN1_VS_TOF_PIN1_I2S_I2N_CORR, pin1_i2pos1_cor_tof, pin1);

        }

        if (hasImplant){
            if (pin1_i2pos1_cor_tof < 1012 && pin1_i2pos1_cor_tof > 917) {
                if (pin1 < 528 && pin1 > 494)
                    has72Co = true;
            }
        }
    }


    //--------------------- PLOTTING PSPMT POSITION --------------------

    if (hasImplant)
        plot(DD_IMPLANT_PSPMT, position.first * pspmtScale + pspmtOffset, position.second * pspmtScale + pspmtOffset);

    if (hasIon)
        plot(DD_ION_PSPMT, position.first * pspmtScale + pspmtOffset, position.second * pspmtScale + pspmtOffset);

    if (hasDecay)
        plot(DD_DECAY_PSPMT, position.first * pspmtScale + pspmtOffset, position.second * pspmtScale + pspmtOffset);


    //------------------- Correlating Decays with Implants ---------------

    bool has72CoDecay = false;
    double dist;

    static double decay_window = 1 * 59 * pow(10.0, -3.0) / Globals::get()->GetEventLengthInSeconds();
    static double pixel_time[Px][Py] = {};
    if (pixel_time[11][11] == 0) {
        Initialize_Array(pixel_time, decay_window);
    }

    for (int n = 0; n < Px; n++){
        for (int m = 0; m < Py; m++){
            if (pixel.first == n && pixel.second == m) {
                if (has72Co){
                    pixel_time[n][m] = 0;
                    if (implant_map.count(pixel) > 0){
                        implant_map.erase(pixel);
                    }
                    implant_map.insert(make_pair(pixel, position));
                }
                else if (hasDecay && pixel_time[n][m] < decay_window) {  //only look for decays in the decay window
                    has72CoDecay = true;
                    plot(D_DECAY_TIME, pixel_time[n][m] * Globals::get()->GetEventLengthInSeconds() * 10000);
                    //pixel_time[n][m] = decay_window;
                }
            }
            pixel_time[n][m]++;
        }
    }
    if (has72CoDecay) {
        plot(DD_CORR_DECAY_PSPMT, position.first * pspmtScale + pspmtOffset, position.second * pspmtScale + pspmtOffset);
        plot(DD_CORR_DECAY_PIXEL, pixel.first + 1, pixel.second + 1);
        dist = pow( pow( implant_map.at(pixel).first - position.first,2) +
                           pow( implant_map.at(pixel).second - position.second, 2), 0.5);
        plot(D_DECAY_DIST, dist * 5000);

    }

    //------------------- Plotting Correlated decays with particular gates ---------------------

            for (vector<ChanEvent *>::const_iterator iterator2 = geEvts.begin();
             iterator2 != geEvts.end(); iterator2++) {

                if (hasDecay) {
                    plot(D_GE_CORR_DECAY, (*iterator2)->GetCalibratedEnergy());
                }
              //  if (!hasVeto) {
                    plot(D_GE_DECAY_GATED, (*iterator2)->GetCalibratedEnergy());
              //  }
                if (has72CoDecay)
                    plot(D_GE_SELECT_GATED, (*iterator2)->GetCalibratedEnergy());

    }

return (true);

}

double E14060Processor::CorrectToFByI2Pos(const std::string &name,
					  const double &tof,
					  const double &i2ns) {

  double slope = 0.0;
  double intercept = 0.0;
  double dammOffset = 1000.;

  if (name == "pin1_i2n") {
    slope = -2.298;
    intercept = 5056;
  } else if (name == "pin2_i2n"){
    slope = -2.298;
    intercept = 5056;
  }
  return tof - slope * i2ns - intercept + dammOffset;
}

void E14060Processor::Initialize_Array(double array[Px][Py], double &val){
    for (int i = 0; i < Px; i++){
        for (int j = 0; j < Py; j++){
            array[i][j] = val;
        }
    }

}