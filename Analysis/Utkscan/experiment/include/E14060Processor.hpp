/** \file E14060Processor.hpp
 * \brief Experiment specific processor to handle the e14060 experiment at
 * the NSCL.
 *\author S. V. Paulauskas
 *\date September 15, 2016
 */
#ifndef __E14060PROCESSOR_HPP_
#define __E14060PROCESSOR_HPP_
#include <map>

#include "EventProcessor.hpp"
#include "PspmtStruct.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

/// Class to analyze data from e14060
class E14060Processor : public EventProcessor {
public:
    ///Constructor that will accept a user defined range of Clover energys to
    /// gate the PID.
    E14060Processor(std::pair<double, double> &energyRange);

    /** Default Destructor */
    ~E14060Processor();

    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);

    /** Process the event
    * \param [in] event : the event to process
    * \return Returns true if the processing was successful */
    virtual bool Process(RawEvent &event);

private:
    std::pair<double, double> energyRange_; ///!< Ge range ge for cuts on PID


    static const int Px = 12; static const int Py = 12;
    double start_time;
    double decay_window;
    double decay_time;
    double pixel_time[2*Px][2*Py];
    static std::vector<PspmtEvent> past_events;
    double trcNum;


    ///@brief A method that will plot the PID spectra given the inputs.
    ///@param[in] map : The map to search for the desired energy
    ///@param[in] key : The key to search for
    ///@return Returns the associated energy of the TAC or PIN, or returns
    /// zero if it cannot be found. 

    double CorrectTofByPos(const double &i2ns,
			     const double &tof);

    ///@brief A method that will set the types associated with this
    /// processing class
    void SetAssociatedTypes();

//protected:
    TFile *rootfile;
    TTree *roottree;
    double  low_xa;
    double  low_xb;
    double  low_ya;
    double  low_yb;
    double  hi_xa;
    double  hi_xb;
    double  hi_ya;
    double  hi_yb;
    std::string eventType;
    double timestamp;
    PspmtEvent current_event;
    std::vector<PspmtEvent> pastEvents;
    std::vector<double> gammaEvents;
};

#endif
