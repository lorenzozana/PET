//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  GAMOS software  is  copyright of the Copyright  Holders  of *
// * the GAMOS Collaboration.  It is provided  under  the  terms  and *
// * conditions of the GAMOS Software License,  included in the  file *
// * LICENSE and available at  http://fismed.ciemat.es/GAMOS/license .*
// * These include a list of copyright holders.                       *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GAMOS collaboration.                       *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the GAMOS Software license.           *
// ********************************************************************
//
#include "GmDoubleBackToBackParticleSource.hh"
#include "CLHEP/Random/RandFlat.h"

#include "GamosCore/GamosGenerator/include/GmGenerVerbosityMgr.hh"
#include "GamosCore/GamosGenerator/include/GmGenerDistEnergyConstant.hh"
#include "GamosCore/GamosBase/Base/include/GmParameterMgr.hh"
#include "GamosCore/GamosGenerator/include/GmGeneratorDistributionFactories.hh"

#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "GmGenerDistTimeConstant.hh"
#include "GmGenerDistDirectionRandom.hh"
#include "GmGenerDistPositionPoint.hh" 
#include "GmGenerDistEnergyConstant.hh"

#include "Randomize.hh"

//-----------------------------------------------------------------------
GmDoubleBackToBackParticleSource::GmDoubleBackToBackParticleSource( const G4String& name, const G4String& partName, const G4double ener ): GmParticleSource( name )
{
  theType = "DoubleBackToBackParticleSource";



  theTimeDistribution = new GmGenerDistTimeConstant;
  theActivity = 1.*becquerel; 

  theDirectionDistribution = new GmGenerDistDirectionRandom;

  thePositionDistribution = new GmGenerDistPositionPoint;  // source centered at (0,0,0)

  //  if( name == "F18" && energyDist == "BetaDecay" ) {
  theEnergy = ener;
  if( bBiasDistributions ) {
    BiasEnergy();
  }
  theEnergyDistribution = new GmGenerDistEnergyConstant;
  static_cast<GmGenerDistEnergyConstant*>(theEnergyDistribution)->SetEnergy(theEnergy);  // could be extracted from source info at Generate(), but it is more efficient here
  CheckDistributionsExist();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  theParticleDef = partTable->FindParticle( partName);
  if( !theParticleDef ) G4Exception("GmDoubleBackToBackParticleSource  " + GetName() + "  particle does not exist: " + partName);
}


//-----------------------------------------------------------------------
G4PrimaryVertex* GmDoubleBackToBackParticleSource::GenerateVertex( G4double time )
{
#ifndef GAMOS_NO_VERBOSE
  if( GenerVerb(infoVerb) ) G4cout << " G4PrimaryVertex* GmDoubleBackToBackParticleSource::GenerateVertex( G4double time ) " << G4endl;
#endif
  thePosition = thePositionDistribution->GeneratePosition( this );
  if( bBiasDistributions ) {
    BiasPosition();
  }

#ifndef GAMOS_NO_VERBOSE
  if( GenerVerb(infoVerb) ) G4cout << "GmDoubleBackToBackParticleSource::GenerateVertex pos " << thePosition << G4endl;
#endif

  G4PrimaryVertex* vtx = new G4PrimaryVertex( thePosition, time );

  theEnergy = theEnergyDistribution->GenerateEnergy( this );
#ifndef GAMOS_NO_VERBOSE
  if( GenerVerb(infoVerb) ) G4cout << "GmDoubleBackToBackParticleSource::GenerateVertex ener " << theEnergy << G4endl;
#endif

  theDirection = theDirectionDistribution->GenerateDirection( this );
  theDirection /= theDirection.mag();
  if( bBiasDistributions ) {
    BiasDirection();
  }

  // //***********************************************
  // // Original Momentum Direction Code
  // //***********************************************
  
  //   G4ThreeVector mom = theDirection * std::sqrt(theEnergy*theEnergy + 2.*theParticleDef->GetPDGMass()*theEnergy);
  // #ifndef GAMOS_NO_VERBOSE
  //   if( GenerVerb(infoVerb) ) G4cout << " GmDoubleBackToBackParticleSource::GenerateVertex  mom " << mom << " theEnergy " << theEnergy << G4endl;
  // #endif
  //   G4PrimaryParticle* particle = new G4PrimaryParticle( theParticleDef, mom.x(), mom.y(), mom.z() );
  //   vtx->SetPrimary(particle);
  
  //   G4PrimaryParticle* particle2 = new G4PrimaryParticle( theParticleDef, -mom.x(), -mom.y(), -mom.z() );
  //   vtx->SetPrimary(particle2);
  
  // //***********************************************
  // // Polarization Code Efforts
  // //***********************************************
  
  // // Adding polarization to a double back to back gamma source 
  // //WARNING Polarization is hardcoded so only use gamma sources

  // // if( partName == "gamma") {
    
  // G4ThreeVector ranv(1,0,0);
  // G4double randphi= -pi+G4UniformRand()*2*pi; //rand phi in range -pi to pi
  // ranv.setPhi(randphi);
    
  // // particleSource->SetCurrentSourceto (0);
  // G4double xPol = ranv.x();
  // G4double yPol = ranv.y();
  // // G4double zPol = ranv.z();

  // particle->SetPolarization(xPol, yPol, 0);
  // // particle->SetPolarization(1,0,0);
  // // particleSource->SetCurrentSourceto (1);
  
  // // G4cout << " ranv " << ranv << G4endl;
  // // G4cout << " randphi " << randphi << G4endl;
  
  // //G4ThreeVector ranv2(-1,0,0);
  // //G4double randphi2= -pi+G4UniformRand()*2*pi; //rand phi in range -pi to pi
  
  // //   ranv2.setPhi(randphi2);
    
  // ranv.rotateZ(90.*deg);
  
  // G4double xPol2 = ranv.x();
  // G4double yPol2 = ranv.y();
  // //  G4double zPol2 = ranv.z();
 
  // //particle2->SetPolarization(-ranv);
  // particle2->SetPolarization(-xPol2, -yPol2, 0);
  // // particle2->SetPolarization(-1,0,0);
  // // // }
  
 
  //  //New Polarization code to deal properly with polarization

  //    G4ThreeVector ranv(1,0,0);

  //  // G4ThreeVector ranv = mom.orthogonal();

  //  // G4cout << "Momemtum " << mom.x() << " "<< mom.y() << " " << mom.z() << G4endl;
  //  // G4cout << "ranv " << ranv.x() << " "<< ranv.y() << " " << ranv.z() << G4endl;


  //  // G4double dotTest = mom.dot(ranv);

  // // G4cout << "Dot Product " << dotTest  << G4endl;

  //  G4ThreeVector ranvUnit = ranv.unit();

  //  G4double randphi= -pi+G4UniformRand()*2*pi; //rand phi in range -pi to pi
  
  //  ranvUnit.setPhi(randphi);

  //  G4double xPol = ranvUnit.x();
  //  G4double yPol = ranvUnit.y();
  
  //  //IMPORTANT!! Particle is going in the negative direction and is hence the second photon in the analysis
  
  //  particle->SetPolarization(xPol, yPol, 0);
  //  // particle->SetPolarization(0,1,0);
  
  //  // ranv.rotateZ(90.*deg);

  //  ranvUnit.rotate(mom, 90.*deg);
  
  //  //IMPORTANT!! Particle2 is going in the positive direction and is hence the second photon in the analysis
  
  //  xPol = ranvUnit.x();
  //  yPol = ranvUnit.y();
  
  //  particle2->SetPolarization(xPol, yPol, 0);
  //  //  particle2->SetPolarization(1,0,0);
  
  //*****************************************************
  // Polarization implimentation Version 3.0!
  //*****************************************************
  // Generate random photon direction

  // G4double cost = 1. - 2.*G4UniformRand();

  // G4double sint = sqrt((1.-cost)*(1.+cost));

  // G4double phi = 2.*Pi*G4UniformRand();

  // G4double sinp = sin(phi);

  // G4double cosp = cos(phi);

  // // Determine polarization of new photon

  // G4double sx = cost*cosp;

  // G4double sy = cost*sinp;

  // G4double sz = -sint;

  // G4ThreeVector photonPolarization(sx, sy, sz);

  // G4ThreeVector perp = mom.cross(photonPolarization);

  // photonPolarization = cosp * photonPolarization + sinp * perp;

  // photonPolarization = photonPolarization.unit();

  // //IMPORTANT!! Particle is going in the negative direction and is hence the second photon in the analysis
  
  //  G4double xPol = photonPolarization.x();
  //  G4double yPol = photonPolarization.y();
  //  G4double zPol = photonPolarization.z();

  //  // particle->SetPolarization(xPol, yPol, zPol);
  //  //particle->SetPolarization(photonPolarization);


  //   particle->SetPolarization(0,1,0);
		       	
  //  //IMPORTANT!! Particle2 is going in the positive direction and is hence the second photon in the analysis

  //  photonPolarization =

  //    mom.cross(photonPolarization).unit();

  //  G4double xPol2 = photonPolarization.x();
  //  G4double yPol2 = photonPolarization.y();
  //  G4double zPol2 = photonPolarization.z();


  //  // particle2->SetPolarization(xPol2, yPol2, zPol2);
  //  // particle2->SetPolarization(photonPolarization);
  //    particle2->SetPolarization(1,0,0);

  
  //*****************************************************
  // Polarization implimentation Version 4.0!
  //*****************************************************

  // Generate random photon direction
  
  G4double cost = 1. - 2.*G4UniformRand();
  
  G4double sint = sqrt((1.-cost)*(1.+cost));
  
  G4double phi = twopi*G4UniformRand();
  
  G4double sinp = sin(phi);
  
  G4double cosp = cos(phi);
  
  G4double px = sint*cosp;
  
  G4double py = sint*sinp;
  
  G4double pz = cost;
  
  // Create photon momentum direction vector
  
  G4ThreeVector photonDirection(px, py, pz);

  G4ThreeVector mom = photonDirection * std::sqrt(theEnergy*theEnergy + 2.*theParticleDef->GetPDGMass()*theEnergy);
#ifndef GAMOS_NO_VERBOSE
  if( GenerVerb(infoVerb) ) G4cout << " GmDoubleBackToBackParticleSource::GenerateVertex  mom " << mom << " theEnergy " << theEnergy << G4endl;
#endif
  G4PrimaryParticle* particle = new G4PrimaryParticle( theParticleDef, mom.x(), mom.y(), mom.z() );
  vtx->SetPrimary(particle);

  G4PrimaryParticle* particle2 = new G4PrimaryParticle( theParticleDef, -mom.x(), -mom.y(), -mom.z() );
  vtx->SetPrimary(particle2);


  // Determine polarization of new photon

  G4double sx = cost*cosp;

  G4double sy = cost*sinp;

  G4double sz = -sint;

  G4ThreeVector photonPolarization(sx, sy, sz);

  G4ThreeVector perp = photonDirection.cross(photonPolarization);

  phi = twopi*G4UniformRand();

  sinp = sin(phi);

  cosp = cos(phi);

  photonPolarization = cosp * photonPolarization + sinp * perp;

  photonPolarization = photonPolarization.unit();


  //IMPORTANT!! Particle is going in the negative direction and is hence the second photon in the analysis
  
  G4double xPol = photonPolarization.x();
  G4double yPol = photonPolarization.y();
  G4double zPol = photonPolarization.z();

  //particle->SetPolarization(xPol, yPol, zPol);
  // particle->SetPolarization(photonPolarization);
  particle->SetPolarization(0,0,0);
  
  //IMPORTANT!! Particle2 is going in the positive direction and is hence the first photon in the analysis
  
  photonPolarization = photonDirection.cross(photonPolarization).unit();
  
  G4double xPol2 = photonPolarization.x();
  G4double yPol2 = photonPolarization.y();
  G4double zPol2 = photonPolarization.z();
  
  
  // particle2->SetPolarization(xPol2, yPol2, zPol2);
  // particle2->SetPolarization(photonPolarization);
  particle2->SetPolarization(0,0,0);
  
  
  // G4cout << " Gamma1 Polarization " << checkX << " "<< checkY << " "<< checkZ <<" Gamma2 Polarization " << checkX2 <<" "<< checkY2 <<" "<< checkZ2 << G4endl;
  
  return vtx;
  
}
