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

/// Class to analyze data from e14060
class E14060Processor : public EventProcessor {
public:
    ///Constructor that will accept a user defined range of Clover energys to
    /// gate the PID.
    E14060Processor(std::pair<double, double> &energyRange);

    /** Default Destructor */
    ~E14060Processor() {};

    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);

    /** Process the event
    * \param [in] event : the event to process
    * \return Returns true if the processing was successful */
    virtual bool Process(RawEvent &event);

private:
    std::pair<double, double> energyRange_; ///!< Ge range ge for cuts on PID

    static const int Px = 24; static const int Py = 24;
    double decay_window;
    double pixel_time[Px][Py];

    ///@brief A method that will plot the PID spectra given the inputs.
    ///@param[in] map : The map to search for the desired energy
    ///@param[in] key : The key to search for
    ///@return Returns the associated energy of the TAC or PIN, or returns
    /// zero if it cannot be found. 
    double FindPinOrTacEnergy(const std::map<std::string,double> &mp, 
		 const std::string &key);

    double CorrectToFByI2Pos(const std::string &name, const double &tof,
			     const double &i2ns);

    ///@brief A method that will set the types associated with this
    /// processing class
    void SetAssociatedTypes();

    void Initialize_Array(double array[Px][Py], double &val);
};

#endif
