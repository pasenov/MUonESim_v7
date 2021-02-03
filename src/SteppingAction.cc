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
// $Id: SteppingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "StackingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4VTrajectory.hh"
#include "G4VProcess.hh"
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
#include <math.h>
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* DA, RunAction* RA, EventAction* EA, StackingAction* SA)
:G4UserSteppingAction(), fDetectorconstruction(DA), fRunaction(RA), fEventaction(EA), fStackingaction(SA)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 //get volume of the current step
  G4VPhysicalVolume* volume 
  = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  //step size
  //G4double stepSize = step->GetStepLength();  

  //collect energy step by step
  //G4double edep = step->GetTotalEnergyDeposit();

  //track position coordinates
  G4StepPoint* pre = step->GetPreStepPoint();
  G4StepPoint* post = step->GetPostStepPoint();

  G4double xp,yp,zp,x,y,z;
  xp = pre->GetPosition().x();
  yp = pre->GetPosition().y();
  zp = pre->GetPosition().z();
  x = post->GetPosition().x();
  y = post->GetPosition().y();
  z = post->GetPosition().z();

  //geometry parameters
  /*G4double TargetSizeX, TargetSizeY, TargetSizeZ;
  TargetSizeX = fDetectorconstruction->GetSizeXT();
  TargetSizeY = fDetectorconstruction->GetSizeYT();
  TargetSizeZ = fDetectorconstruction->GetSizeZT();*/

  //momentum direction of primary when exiting target
  G4double momDirPx, momDirPy, momDirPz;

  //momentum direction of secondary when exiting target
  G4double momDirSx, momDirSy, momDirSz;

  //kinetic energies when exiting target
  G4double ekinPre;
  G4double ekinPost;
  
  //charge of secondaries
  G4double charge;

  //get track length, track ID, track length, global time and track's vertex position of the current step
  G4Track* track = step->GetTrack();
  /*G4double length = track->GetTrackLength();
  G4double globalTime = track->GetGlobalTime();*/
  G4ThreeVector primVert = track->GetVertexPosition();
  //G4int trackID = track->GetTrackID();
  
  //process names and ID
  G4String processName;
  G4int processID = 0;

  //exit points after target for primary and secondaries
  G4ThreeVector pointExPrimT, pointExSecT;

  //exit angles of primary and secondaries from target
  G4double anglePrim, angleSec;

  //particle definition
  const G4ParticleDefinition* particleDefinition = track->GetParticleDefinition();
  
  //analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);

  //track length of primary particle
  if (track->GetTrackID() == 1)  {

    ekinPre = pre->GetKineticEnergy();
    ekinPost = post->GetKineticEnergy();
    /*G4cout
       << "\n ekinPre = "
       << G4BestUnit(ekinPre, "Energy")
       << "\n ekinPost = "
       << G4BestUnit(ekinPost, "Energy")
       << "\n z = "
       << G4BestUnit(z, "Length")
       << G4endl;*/
    
    //kinetic energy of primary particle at the exit of the target
    if (volume == fDetectorconstruction->GetTarget())   {

	if (post->GetStepStatus() == fGeomBoundary)   {

	   ekinPost = post->GetKineticEnergy();
	   pointExPrimT = post->GetPosition();

	   momDirPx = track->GetMomentumDirection().x();
	   momDirPy = track->GetMomentumDirection().y();
	   momDirPz = track->GetMomentumDirection().z();

 	   anglePrim = acos(0*momDirPx + 0*momDirPy + 1*momDirPz);

           /*G4cout 
              << "\n Primary has exited target at: " 
              << G4BestUnit(pointExPrimT, "Length")
              << "\n with momentum direction: " 
              << momDirPx
              << " "
              << momDirPy
              << " "
              << momDirPz
	      << "\n and energy: "
	      << G4BestUnit(ekinPost, "Energy")
	      << "\n and angle: "
	      << anglePrim
              << G4endl; */

	   if (anglePrim < 30*mrad)   {
	      fEventaction->AddAnglePrimary(anglePrim);
	      fEventaction->AddEnergyExitTargetPrimary(ekinPost);
	   }

        }

    }

  }

  //kinetic energy of secondaries at the exit of the target
  if (track->GetParentID() == 1)  {

	if (post->GetStepStatus() == fGeomBoundary)   {

	   ekinPost = post->GetKineticEnergy();
	   pointExSecT = post->GetPosition();
	   processName = track->GetCreatorProcess()->G4VProcess::GetProcessName();
	   
	   momDirSx = track->GetMomentumDirection().x();
	   momDirSy = track->GetMomentumDirection().y();
	   momDirSz = track->GetMomentumDirection().z();

 	   angleSec = acos(0*momDirSx + 0*momDirSy + 1*momDirSz);
	   
	   charge = post->GetCharge();

	   if (track->GetLogicalVolumeAtVertex() == fDetectorconstruction->GetLTarget())  {
	   
	      if(charge != 0)   {
	         if ((angleSec < 30*mrad) && (ekinPost > 0.0*GeV) ) {
	            fEventaction->AddAngleCharged(angleSec);
	      	    fEventaction->AddEnergyExitTargetCharged(ekinPost);
	      	    fEventaction->AddNbTracks();
	      	    
	      	    if (processName == "muIoni")  {
	      	       fEventaction->AddAngleIoni(angleSec);
	      	       fEventaction->AddEnergyExitTargetIoni(ekinPost);
	      	    }
	      	    
	      	    if (processName == "muPairProd")  {
	      	       fEventaction->AddAnglePP(angleSec);
	      	       fEventaction->AddEnergyExitTargetPP(ekinPost);
	      	    }
	      	 }
	         if ((angleSec < 30*mrad) && (ekinPost > 0.5*GeV) ) {
	            fEventaction->AddAngleChargedCut05GeV(angleSec);
	      	    fEventaction->AddEnergyExitTargetChargedCut05GeV(ekinPost);
	      	    fEventaction->AddNbTracksCut05GeV();
	      	    
	      	    if (processName == "muIoni")  {
	      	       fEventaction->AddAngleIoniCut05GeV(angleSec);
	      	       fEventaction->AddEnergyExitTargetIoniCut05GeV(ekinPost);
	      	    }
	      	    
	      	    if (processName == "muPairProd")  {
	      	       fEventaction->AddAnglePPCut05GeV(angleSec);
	      	       fEventaction->AddEnergyExitTargetPPCut05GeV(ekinPost);
	      	    }
	      	 }
	      	 if ((angleSec < 30*mrad) && (ekinPost > 1.0*GeV) ) {
	      	    fEventaction->AddAngleChargedCut1GeV(angleSec);
	      	    fEventaction->AddEnergyExitTargetChargedCut1GeV(ekinPost);
	      	    fEventaction->AddNbTracksCut1GeV();
	      	    
	      	    if (processName == "muIoni")  {
	      	       fEventaction->AddAngleIoniCut1GeV(angleSec);
	      	       fEventaction->AddEnergyExitTargetIoniCut1GeV(ekinPost);
	      	    }
	      	    
	      	    if (processName == "muPairProd")  {
	      	       fEventaction->AddAnglePPCut1GeV(angleSec);
	      	       fEventaction->AddEnergyExitTargetPPCut1GeV(ekinPost);
	      	    }
	      	 }
	         
	      }

	      if(particleDefinition == G4Positron::Definition())   {
                 /*G4cout 
                    << "\n Secondary positron has exited target at: " 
                    << G4BestUnit(pointExSecT, "Length")
                    << "\n with momentum direction: " 
                    << momDirSx
                    << " "
                    << momDirSy
                    << " "
                    << momDirSz
                    << G4endl; */

	         if ((angleSec < 30*mrad) && (ekinPost > 1.0*GeV) ) {
	            fEventaction->AddAngleSecondaryPositron(angleSec);
	            fEventaction->AddEnergyExitTargetSecondaryPositron(ekinPost);
	            
	            if (processName == "muIoni") processID = 1; //ionization
	            if (processName == "muPairProd") processID = 2;//pair production
	            if (processName == "muBrems") processID = 3; //Bremsstrahlung
	            if (processName == "muNucl") processID = 4; //nuclear interaction
	            analysisManager->FillH1(11, processID);
	            analysisManager->FillNtupleDColumn(11, processID);
	            analysisManager->AddNtupleRow();
	            
	            /*G4cout
                       << "\n Positron created due to: "
                       << processName
                       << G4endl;*/
	         }

	      }

	      if(particleDefinition == G4Electron::Definition())   {
                 /*G4cout 
                    << "\n Secondary electron has exited target at: " 
                    << G4BestUnit(pointExSecT, "Length")
                    << "\n with momentum direction: " 
                    << momDirSx
                    << " "
                    << momDirSy
                    << " "
                    << momDirSz
		    << "\n and energy: "
		    << G4BestUnit(ekinPost, "Energy")
	      	    << "\n and angle: "
	            << angleSec
                    << G4endl; */

	         if ((angleSec < 30*mrad) && (ekinPost > 1.0*GeV) ) {
	            fEventaction->AddAngleSecondaryElectron(angleSec);
	            fEventaction->AddEnergyExitTargetSecondaryElectron(ekinPost);
	            
	            if (processName == "muIoni") processID = 1; //ionization
	            if (processName == "muPairProd") processID = 2;//pair production
	            if (processName == "muBrems") processID = 3; //Bremsstrahlung
	            if (processName == "muNucl") processID = 4; //nuclear interaction
	            analysisManager->FillH1(11, processID);
	           
	            /*G4cout
                       << "\n Electron created due to: "
                       << processName
                       << G4endl;*/
	         }

	      }

	      if(particleDefinition == G4MuonPlus::Definition())   {
                 /*G4cout 
                    << "\n Secondary muon plus has exited target at: " 
                    << G4BestUnit(pointExSecT, "Length")
                    << "\n with momentum direction: " 
                    << momDirSx
                    << " "
                    << momDirSy
                    << " "
                    << momDirSz
                    << G4endl; */

	         if ((angleSec < 30*mrad) && (ekinPost > 1.0*GeV) ) {
	            fEventaction->AddAngleSecondaryMuonPlus(angleSec);
	            fEventaction->AddEnergyExitTargetSecondaryMuonPlus(ekinPost);
	         }

	      }

	      if(particleDefinition == G4MuonMinus::Definition())   {
                 /*G4cout 
                    << "\n Secondary muon minus has exited target at: " 
                    << G4BestUnit(pointExSecT, "Length")
                    << "\n with momentum direction: " 
                    << momDirSx
                    << " "
                    << momDirSy
                    << " "
                    << momDirSz
                    << G4endl; */

	         if ((angleSec < 30*mrad) && (ekinPost > 1.0*GeV) ) {
	            fEventaction->AddAngleSecondaryMuonMinus(angleSec);
	            fEventaction->AddEnergyExitTargetSecondaryMuonMinus(ekinPost);
	         }

	      }

	   }

	}

  }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

