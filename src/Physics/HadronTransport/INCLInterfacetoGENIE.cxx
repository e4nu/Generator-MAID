#ifdef GOPT_ENABLE_INCL

//#include "G4INCLConfig.hh"
//#include "G4AblaInterface.hh"

//#include "INCLInterfacetoGENIE.hh"
#include "INCLInterfacetoGENIE.h"

//#include "G4INCLParticleSpecies.hh"
//#include "G4INCLConfigEnums.hh"
//#include <vector>
//#include "G4INCLConfig.hh"

#include "INCLConfigParser.hh"

//#include <iostream>
//#include <iomanip>
//#include <string>
//#include <fstream>
//#include <ostream>
//#include <sstream>
//#include <cassert>
//#include <cstdlib>
//#include <map>
//#include <cstdlib>
//#include <sstream>

using namespace G4INCL;
G4ThreadLocal INCLInterfacetoGENIE *INCLInterfacetoGENIE::theInstance = NULL;

INCLInterfacetoGENIE::INCLInterfacetoGENIE():theINCLModel(NULL) {
}

INCLInterfacetoGENIE::~INCLInterfacetoGENIE() {
  delete theINCLModel;
}

INCLInterfacetoGENIE *INCLInterfacetoGENIE::GetInstance() {
  if ( ! theInstance )
    theInstance = new INCLInterfacetoGENIE;
  return theInstance;
}

void INCLInterfacetoGENIE::DeleteModel() {
  delete theINCLModel; theINCLModel=NULL;
}

G4INCL::Config *INCLInterfacetoGENIE::setINCLConfig() {
  DeleteModel(); // in case the Config is modified

  // can't legally initialize "char *" from string constant
  const char * test[] = {"NULL","-pp","-tFe56","-N100","-E100","-dabla07"};
  int argcc = 6;
  INCLConfigParser theParser;

  // RWH throw away above const-ness
  G4INCL::Config *aConfig=theParser.parse(argcc,(char**)test);
  std::cerr << "==RWH== INCLInterfacetoGENIE create theConfig" << std::endl;

  theConfig = aConfig;
  return theConfig;

}

G4INCL::INCL *INCLInterfacetoGENIE::GetINCLModel() {
  if ( ! theINCLModel ) {

    //   G4INCL::Config *aConfig = new G4INCL::Config(theConfig);
    std::cerr << "==RWH== INCLInterfacetoGENIE::GetINCLModel new G4INCL::INCL" << std::endl;
    setINCLConfig();
    theINCLModel = new G4INCL::INCL(theConfig);
    // ownership of the aConfig object is taken over by the INCL model engine
  }
  return theINCLModel;
}

#endif // __GENIE_INCL_ENABLED__
