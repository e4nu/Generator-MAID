//____________________________________________________________________________
/*
  Copyright (c) 2023-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
  University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Resonance/XSection/MAIDRESVectFormFactorsEMn.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDRESVectFormFactorsEMn::MAIDRESVectFormFactorsEMn() :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMn")
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMn::MAIDRESVectFormFactorsEMn(string config) :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMn", config)
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMn::~MAIDRESVectFormFactorsEMn()
{

}
//____________________________________________________________________________
RESVectFFAmplitude MAIDRESVectFormFactorsEMn::Compute( const Interaction interaction ) {
  RESVectFFAmplitude ampl ; 

  const InitialState & init_state = interaction.InitState();
  const Target & target = init_state.Tgt();
  const Resonance_t res = interaction.ExclTag().Resonance();
  
  //Get kinematical parameters
  const Kinematics & kinematics = interaction.Kine();
  double q2 = kinematics.q2();

  if( res == kP33_1232 ) { 
    double Mnuc2  = TMath::Power(target.HitNucMass(), 2);
    double Mnuc3  = TMath::Power(target.HitNucMass(), 3);
    double MR  = utils::res::Mass(res);
    double MR2 = TMath::Power(MR,2);
    double W = kinematics.W();
    double W2 = TMath::Power(W,2);
    double k = TMath::Sqrt( TMath::Power((W2-Mnuc2+q2)*0.5/W,2) -q2 );
    double kw = (W2-Mnuc2)*0.5/W ; // k(W,0);
    double kR = (MR2-Mnuc2)*0.5/MR; // kR=k(W=MR,0)
    double kR2 = TMath::Power(kR,2);
    double kappaR = (MR2 - Mnuc2)*0.5/MR;                                 // kR(W=MR,Q=0)
    double cdelta = TMath::Sqrt(Mnuc3*kappaR/(4*kPi*kAem*MR*kR2));

    // Dipole form facor
    double GD = 1./TMath::Power(1-q2/fDipoleMass,2) ;
    // Transition form factors
    double AM = fAM0_P33_1232 * (1. - fBetaM_P33_1232 * q2) * TMath::Exp( fGammaM_P33_1232 * q2 ) * GD / 1000.;//* TMath::Sqrt(115./130.) ;* (k/kw) *
    double AE = fAE0_P33_1232 * (1. - fBetaE_P33_1232 * q2) * TMath::Exp( fGammaE_P33_1232 * q2 ) * GD / 1000.;//* TMath::Sqrt(115./130.) ;* (k/kw) *
    double AS = fAC0_P33_1232 * (1. - fBetaC_P33_1232 * q2) / (1. - fDC_P33_1232 * q2/(4.*Mnuc2))*TMath::Exp( fGammaC_P33_1232 * q2 ) * GD / 1000. ; 
    AS *= TMath::Power(k,2) / (kw * kR);
    ampl.SetAmplA12( - 0.5 *(3*AE+AM) );
    ampl.SetAmplA32( TMath::Sqrt(3.)*0.5*(AE-AM) );
    ampl.SetAmplS12( - TMath::Sqrt(2.) * AS );
  } else {
    double A120 = fA120N[res] ;
    double A12Alpha = fA12AlphaN[res] ;
    double A12Beta = fA12BetaN[res] ;
    double A320 = fA320N[res] ;
    double A32Alpha = fA32AlphaN[res] ;
    double A32Beta = fA32BetaN[res] ;
    double S120 = fS120N[res] ;
    double S12Alpha = fS12AlphaN[res] ;
    double S12Beta = fS12BetaN[res] ;

    ampl.SetAmplA12( A120 * ( 1 - A12Alpha * q2 ) * TMath::Exp( A12Beta * q2 ) ) ; 
    ampl.SetAmplA32( A320 * ( 1 - A32Alpha * q2 ) * TMath::Exp( A32Beta * q2 ) ) ; 
    ampl.SetAmplS12( S120 * ( 1 - S12Alpha * q2 ) * TMath::Exp( S12Beta * q2 ) ) ; 
  }

  return ampl;
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMn::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMn::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMn::LoadConfig(void)
{
  bool good_config = true ; 

  GetParam( "AM@P33(1232)", fAM0_P33_1232 ) ; 
  GetParam( "AE@P33(1232)", fAE0_P33_1232 ) ;
  GetParam( "AC@P33(1232)", fAC0_P33_1232 ) ; 
  
  GetParam( "BetaM@P33(1232)", fBetaM_P33_1232 ) ;
  GetParam( "BetaE@P33(1232)", fBetaE_P33_1232 ) ;
  GetParam( "BetaC@P33(1232)", fBetaC_P33_1232 ) ;

  GetParam( "GammaM0@P33(1232)", fGammaM_P33_1232 ) ;
  GetParam( "GammaE0@P33(1232)", fGammaE_P33_1232 ) ;
  GetParam( "GammaC0@P33(1232)", fGammaC_P33_1232 ) ;

  GetParam( "DeltaC@P33(1232)", fDC_P33_1232 ) ;
  GetParam( "DipoleMass", fDipoleMass ) ;

  auto kres_list_A12_0_n = GetConfig().FindKeys("A120N@") ;
  if( kres_list_A12_0_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_0_n.begin(); kiter != kres_list_A12_0_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA120N[res_id] ) ; 
  }

  auto kres_list_A12_alpha_n = GetConfig().FindKeys("A12AlphaN@") ;
  if( kres_list_A12_alpha_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_alpha_n.begin(); kiter != kres_list_A12_alpha_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA12AlphaN[res_id] ) ; 
  }

  auto kres_list_A12_beta_n = GetConfig().FindKeys("A12BetaN@") ;
  if( kres_list_A12_beta_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_beta_n.begin(); kiter != kres_list_A12_beta_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA12BetaN[res_id] ) ; 
  }


  auto kres_list_A32_0_n = GetConfig().FindKeys("A320N@") ;
  if( kres_list_A32_0_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_0_n.begin(); kiter != kres_list_A32_0_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA320N[res_id] ) ; 
  }

  auto kres_list_A32_alpha_n = GetConfig().FindKeys("A32AlphaN@") ;
  if( kres_list_A32_alpha_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_alpha_n.begin(); kiter != kres_list_A32_alpha_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA32AlphaN[res_id] ) ; 
  }

  auto kres_list_A32_beta_n = GetConfig().FindKeys("A32BetaN@") ;
  if( kres_list_A32_beta_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_beta_n.begin(); kiter != kres_list_A32_beta_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA32BetaN[res_id] ) ; 
  }


  auto kres_list_S12_0_n = GetConfig().FindKeys("S120N@") ;
  if( kres_list_S12_0_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_0_n.begin(); kiter != kres_list_S12_0_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fS120N[res_id] ) ; 
  }

  auto kres_list_S12_alpha_n = GetConfig().FindKeys("S12AlphaN@") ;
  if( kres_list_S12_alpha_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_alpha_n.begin(); kiter != kres_list_S12_alpha_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fS12AlphaN[res_id] ) ; 
  }

  auto kres_list_S12_beta_n = GetConfig().FindKeys("S12BetaN@") ;
  if( kres_list_S12_beta_n.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_beta_n.begin(); kiter != kres_list_S12_beta_n.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fS12BetaN[res_id] ) ; 
  }

  if( ! good_config ) { 
    LOG("MAIDRESVectFormFactorsEMn", pERROR ) << " Configuration failed.";
    exit(78) ;
  }

}
