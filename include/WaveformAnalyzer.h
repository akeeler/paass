/** \file WaveformProcessor.h
 * 
 * Class for handling Scintillator traces. 
 */

#ifndef __WAVEFORMANALYZER_H_
#define __WAVEFORMANALYZER_H_

#include "Trace.h"
#include "TraceAnalyzer.h"

class WaveformAnalyzer : public TraceAnalyzer
{
 public:
    struct FitData {
	size_t n;
	double * y;
	double * sigma;
	float WID;
	float DKAY;
    };

    WaveformAnalyzer(); // no virtual c'tors
    virtual void DeclarePlots(void);
    virtual void Analyze(Trace &, const std::string &,
			 const std::string &);
    
    virtual ~WaveformAnalyzer() {/* do nothing */};
};

#endif // __WAVEFORMANALYZER_H_
