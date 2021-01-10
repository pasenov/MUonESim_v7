//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm18/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "RunAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTrajectory.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(DetectorConstruction* DA, RunAction* RA, EventAction* EA)
:G4UserStackingAction(), fDetectorconstruction(DA), fRunaction(RA), fEventaction(EA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  //define speed unit
  new G4UnitDefinition ("meter/second", "m/s", "Speed", m/s);
  
  //energy spectrum, polarization, initial momentum direction, track weight, velocity and local time of secondaries and tertiaries
  //
  energy = track->GetKineticEnergy();
  xpolarization = track->GetPolarization().x();
  ypolarization = track->GetPolarization().y();
  zpolarization = track->GetPolarization().z();
  charged = (track->GetDefinition()->GetPDGCharge() != 0.);
  G4VPhysicalVolume* volume = track->GetVolume();
  /*G4double trackWeight = track->GetWeight();
  G4double initVelocity = track->GetVelocity();
  G4double kinEnergy = track->GetKineticEnergy();
  G4double time = track->GetLocalTime();*/

  //geometry parameters
  /*G4double TargetSizeX, TargetSizeY, TargetSizeZ;
  G4double pos;
 
  TargetSizeX=fDetectorconstruction->GetSizeXT();
  TargetSizeY=fDetectorconstruction->GetSizeYT();
  TargetSizeZ=fDetectorconstruction->GetSizeZT();*/

  //G4ThreeVector primVert = track->GetVertexPosition();
  /*G4double x_prim, y_prim, z_prim;
  x_prim = track->GetVertexPosition().x();
  y_prim = track->GetVertexPosition().y();
  z_prim = track->GetVertexPosition().z();*/

  //secret number
  //G4int iSecret;

  //G4cout 
      //<< "\n Local time: " 
      //<< G4BestUnit(time, "Time")
      //<< "\n Initial velocity: " 
      //<< G4BestUnit(initVelocity, "Speed")
      //<< "\n Kinetic energy: "
      //<< G4BestUnit(kinEnergy, "Energy")
      //<< "\n Track weight: " 
      //<< trackWeight
      //<< G4endl;

  //particle definition
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();

  G4ThreeVector initMom = track->GetMomentum();

  //kinetic energy of secondaries calculation 
  if (track->GetParentID() == 1)  {

    if (volume == fDetectorconstruction->GetTarget()) {

      if(particleDefinition == G4Positron::Definition())   {
        fRunaction->AddNbPositronsTarget();
      }

      if(particleDefinition == G4Electron::Definition())   {
        fRunaction->AddNbElectronsTarget();
      }

      if(particleDefinition == G4MuonPlus::Definition())   {
        fRunaction->AddNbMuonsPlusTarget();
      }

      if(particleDefinition == G4MuonMinus::Definition())   {
        fRunaction->AddNbMuonsMinusTarget();
      }

    }

  }

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
