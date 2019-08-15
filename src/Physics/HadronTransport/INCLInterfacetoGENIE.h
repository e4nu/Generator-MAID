#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef INCLInterfacetoGENIE_hh
#define INCLInterfacetoGENIE_hh 1

//#include "G4INCLParticleSpecies.hh"
//#include "G4INCLConfigEnums.hh"
//#include "G4INCLRandomSeedVector.hh"
//#include <iostream>
//#include <string>
//#include <sstream>
//#include <cassert>
#include "G4INCLCascade.hh"
#include "G4INCLVersion.hh"
#include "G4INCLConfig.hh"
//#include <list>
//#include <sstream>

namespace G4INCL {

  class INCLInterfacetoGENIE {
  public:
    static INCLInterfacetoGENIE *GetInstance();
    static G4ThreadLocal INCLInterfacetoGENIE *theInstance;
    INCLInterfacetoGENIE();
    ~INCLInterfacetoGENIE();

    G4INCL::Config *setINCLConfig();
    G4INCL::INCL   *GetINCLModel();

    //private:
    void DeleteModel();
    G4INCL::Config *theConfig;
    G4INCL::INCL *theINCLModel;
  };
}

#endif
#endif // GOPT_ENABLE_INCL
//
