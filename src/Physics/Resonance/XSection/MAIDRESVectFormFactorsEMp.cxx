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
#include "Physics/Resonance/XSection/MAIDRESVectFormFactorsEMp.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDRESVectFormFactorsEMp::MAIDRESVectFormFactorsEMp() :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMp")
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMp::MAIDRESVectFormFactorsEMp(string config) :
  RESVectFormFactorsI("genie::MAIDRESVectFormFactorsEMp", config)
{

}
//____________________________________________________________________________
MAIDRESVectFormFactorsEMp::~MAIDRESVectFormFactorsEMp()
{

}
//____________________________________________________________________________
RESVectFFAmplitude MAIDRESVectFormFactorsEMp::Compute( const Interaction interaction ) {
  RESVectFFAmplitude ampl ; 

  const InitialState & init_state = interaction.InitState();
  const ProcessInfo &  proc_info  = interaction.ProcInfo();
  const Target & target = init_state.Tgt();
  const Resonance_t res = interaction.ExclTag().Resonance();
  
  //Get kinematical parameters
  const Kinematics & kinematics = interaction.Kine();
  double q2 = kinematics.q2();

  if( res == kP33_1232 ) { 
    //Compute auxiliary & kinematical factors
    double W2     = TMath::Power( kinematics.W(), 2);
    double Mnuc2  = TMath::Power(target.HitNucMass(), 2);
    double tau = -q2/(4.*Mnuc2);  
    double MR  = utils::res::Mass(res);
    double kR = (MR*MR - Mnuc2)/(2.*MR);  //LAR choice, k(W=MR,Q2=0) Q2=0 real photon
    double kgcm0 = (W2 - Mnuc2)/(2.*MR); // kW at equation 5
    double egcm = (W2+q2-Mnuc2)/(2.*MR); //photon energy at the W center of mass frame
    double qcm = TMath::Sqrt(pow(egcm,2)-q2); //photon momentum k in Equation 3
    
    double Fq = 1./TMath::Power(1-q2/0.71,2)*qcm/kgcm0;
    double AM = fAM0_P33_1232 * (1. - fBetaM_P33_1232 * q2) * TMath::Exp( fGammaM_P33_1232 * q2 ) * Fq;
    double AE = fAE0_P33_1232 * (1. - fBetaE_P33_1232 * q2) * TMath::Exp( fGammaE_P33_1232 * q2 ) * Fq;
    double AS = fAS0_P33_1232 * (1. - fBetaS_P33_1232 * q2) / (1. + 4.9*tau)*qcm/kR*TMath::Exp( fGammaS_P33_1232 * q2 ) * Fq; // ??
    ampl.SetAmplA12( (3.*AE+AM)/2./1000. );
    ampl.SetAmplA32( TMath::Sqrt(3.)/2.*(AE-AM)/1000. );
    ampl.SetAmplS12( TMath::Sqrt(2.)*AS/1000. );
  } else if ( res == kP11_1440) {
    ampl.SetAmplA12( fA120P[res] * ( 1 + fA12AlphaP[res] * q2 + fA12DeltaP * pow( q2, 4 ) ) * TMath::Exp( fA12BetaP[res] * q2 ) / 1000. );
    ampl.SetAmplS12( fS120P[res] * ( 1 + fS12AlphaP[res] * q2 + fS12DeltaP * pow( q2, 4 ) ) * TMath::Exp( fS12BetaP[res] * q2 ) / 1000. );
    ampl.SetAmplA32( 0. ) ;
  } else { 
    double A120 = fA120P[res] ;
    double A12Alpha = fA12AlphaP[res] ;
    double A12Beta = fA12BetaP[res] ;
    double A320 = fA320P[res] ;
    double A32Alpha = fA32AlphaP[res] ;
    double A32Beta = fA32BetaP[res] ;
    double S120 = fS120P[res] ;
    double S12Alpha = fS12AlphaP[res] ;
    double S12Beta = fS12BetaP[res] ;

    ampl.SetAmplA12( A120 * ( 1 - A12Alpha * q2 ) * TMath::Exp( A12Beta * q2 ) / 1000. ) ; 
    ampl.SetAmplA32( A320 * ( 1 - A32Alpha * q2 ) * TMath::Exp( A32Beta * q2 ) / 1000. ) ; 
    ampl.SetAmplS12( S120 * ( 1 - S12Alpha * q2 ) * TMath::Exp( S12Beta * q2 ) / 1000. ) ; 
  }

  return ampl;
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMp::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMp::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectFormFactorsEMp::LoadConfig(void)
{
  bool good_config = true ; 

  GetParam( "AM0@P33(1232)", fAM0_P33_1232 ) ; 
  GetParam( "AE0@P33(1232)", fAE0_P33_1232 ) ;
  GetParam( "AS0@P33(1232)", fAS0_P33_1232 ) ; 
  
  GetParam( "BetaM@P33(1232)", fBetaM_P33_1232 ) ;
  GetParam( "BetaE@P33(1232)", fBetaE_P33_1232 ) ;
  GetParam( "BetaS@P33(1232)", fBetaS_P33_1232 ) ;

  GetParam( "GammaM@P33(1232)", fGammaM_P33_1232 ) ;
  GetParam( "GammaE@P33(1232)", fGammaE_P33_1232 ) ;
  GetParam( "GammaS@P33(1232)", fGammaS_P33_1232 ) ;

  //  GetParam( "NM@P33(1232)", fNM_P33_1232 ) ;
  //GetParam( "NE@P33(1232)", fNE_P33_1232 ) ;

  GetParam( "A12DeltaP@P11(1440)", fA12DeltaP ) ; 
  GetParam( "S12DeltaP@P11(1440)", fS12DeltaP ) ; 

  auto kres_list_A12_0_p = GetConfig().FindKeys("A120P@") ;
  if( kres_list_A12_0_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_0_p.begin(); kiter != kres_list_A12_0_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA120P[res_id] ) ; 
  }

  auto kres_list_A12_alpha_p = GetConfig().FindKeys("A12AlphaP@") ;
  if( kres_list_A12_alpha_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_alpha_p.begin(); kiter != kres_list_A12_alpha_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA12AlphaP[res_id] ) ; 
  }

  auto kres_list_A12_beta_p = GetConfig().FindKeys("A12BetaP@") ;
  if( kres_list_A12_beta_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A12_beta_p.begin(); kiter != kres_list_A12_beta_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA12BetaP[res_id] ) ; 
  }


  auto kres_list_A32_0_p = GetConfig().FindKeys("A320P@") ;
  if( kres_list_A32_0_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_0_p.begin(); kiter != kres_list_A32_0_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA320P[res_id] ) ; 
  }

  auto kres_list_A32_alpha_p = GetConfig().FindKeys("A32AlphaP@") ;
  if( kres_list_A32_alpha_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_alpha_p.begin(); kiter != kres_list_A32_alpha_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA32AlphaP[res_id] ) ; 
  }

  auto kres_list_A32_beta_p = GetConfig().FindKeys("A32BetaP@") ;
  if( kres_list_A32_beta_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_A32_beta_p.begin(); kiter != kres_list_A32_beta_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fA32BetaP[res_id] ) ; 
  }


  auto kres_list_S12_0_p = GetConfig().FindKeys("S120P@") ;
  if( kres_list_S12_0_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_0_p.begin(); kiter != kres_list_S12_0_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fS120P[res_id] ) ; 
  }

  auto kres_list_S12_alpha_p = GetConfig().FindKeys("S12AlphaP@") ;
  if( kres_list_S12_alpha_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_alpha_p.begin(); kiter != kres_list_S12_alpha_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fS12AlphaP[res_id] ) ; 
  }

  auto kres_list_S12_beta_p = GetConfig().FindKeys("S12BetaP@") ;
  if( kres_list_S12_beta_p.size() == 0 ) good_config = false ; 
  for( auto kiter = kres_list_S12_beta_p.begin(); kiter != kres_list_S12_beta_p.end(); ++kiter ) {
    const RgKey & key = *kiter ;
    vector<string> kv = genie::utils::str::Split(key,"@");
    assert(kv.size()==2);
    Resonance_t res_id = utils::res::FromString( (kv[1]).c_str() );
    GetParam( key, fS12BetaP[res_id] ) ; 
  }

  if( ! good_config ) { 
    LOG("MAIDRESVectFormFactorsEMp", pERROR ) << " Configuration failed.";
    exit(78) ;
  }

}
