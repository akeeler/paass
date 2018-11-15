///@file LinkDef.hpp
///@brief LinkDef for Various ROOT output structures
///@authors T.T. King
///@date 2/6/2018

#ifndef PAASS_LINKDEF_HPP
#define PAASS_LINKDEF_HPP
#ifdef __CINT__
#pragma link C++ struct PspmtEvent+;
#pragma link C++ class std::vector<PspmtEvent>+;
#pragma link C++ struct PidEvent+;
#pragma link C++ class std::vector<PidEvent>+;


#endif

#endif //PAASS_LINKDEF_HPP
