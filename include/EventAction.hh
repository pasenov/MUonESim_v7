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
/// \file electromagnetic/TestEm18/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*, RunAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void   EndOfEventAction(const G4Event*);

    void AddEnergyExitTargetPrimary(G4double edep)            {fEnergyExitTargetPrimary  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryPositron(G4double edep)  {fEnergyExitTargetSecondaryPositron  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryElectron(G4double edep)  {fEnergyExitTargetSecondaryElectron  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryMuonPlus(G4double edep)  {fEnergyExitTargetSecondaryMuonPlus  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetSecondaryMuonMinus(G4double edep) {fEnergyExitTargetSecondaryMuonMinus  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetCharged(G4double edep) {fEnergyExitTargetCharged  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetChargedCut05GeV(G4double edep) {fEnergyExitTargetChargedCut05GeV  += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetChargedCut1GeV(G4double edep) {fEnergyExitTargetChargedCut1GeV  += edep;}; // If angle < 30 mrad
    
    void AddEnergyExitTargetIoni(G4double edep)   {fEnergyExitTargetIoni += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetPP(G4double edep)   {
       if (fNbTracks == 1)  {fEnergyExitTargetPP1 += edep;}
       if (fNbTracks == 2)  {fEnergyExitTargetPP2 += edep;} 
    }; // If angle < 30 mrad
    
    void AddEnergyExitTargetIoniCut05GeV(G4double edep)   {fEnergyExitTargetIoniCut05GeV += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetPPCut05GeV(G4double edep)   {
       if (fNbTracksCut05GeV == 1)  {fEnergyExitTargetPP1Cut05GeV += edep;}
       if (fNbTracksCut05GeV == 2)  {fEnergyExitTargetPP2Cut05GeV += edep;} 
    }; // If angle < 30 mrad
    
    void AddEnergyExitTargetIoniCut1GeV(G4double edep)   {fEnergyExitTargetIoniCut1GeV += edep;}; // If angle < 30 mrad
    void AddEnergyExitTargetPPCut1GeV(G4double edep)   {
       if (fNbTracksCut1GeV == 1)  {fEnergyExitTargetPP1Cut1GeV += edep;}
       if (fNbTracksCut1GeV == 2)  {fEnergyExitTargetPP2Cut1GeV += edep;} 
    }; // If angle < 30 mrad


    void AddAnglePrimary(G4double theta)            	      {if (fBoolPrim == 0)   { fAnglePrimary  += theta; fBoolPrim = 1; } }; // If angle < 30 mrad
    void AddAngleSecondaryPositron(G4double theta)            {if (fBoolSecPositron == 0)   {fAngleSecondaryPositron  += theta; fBoolSecPositron = 1;} }; // If angle < 30 mrad
    void AddAngleSecondaryElectron(G4double theta)            {if (fBoolSecElectron == 0)   {fAngleSecondaryElectron  += theta; fBoolSecElectron = 1;} }; // If angle < 30 mrad
    void AddAngleSecondaryMuonPlus(G4double theta)            {if (fBoolSecMuonPlus == 0)   {fAngleSecondaryMuonPlus  += theta; fBoolSecMuonPlus = 1;} }; // If angle < 30 mrad
    void AddAngleSecondaryMuonMinus(G4double theta)           {if (fBoolSecMuonMinus == 0)  {fAngleSecondaryMuonMinus += theta; fBoolSecMuonMinus = 1;} }; // If angle < 30 mrad
    void AddAngleCharged(G4double theta)            	      {if (fBoolCharged == 0)   { fAngleCharged  += theta; fBoolCharged = 1; } }; // If angle < 30 mrad
    void AddAngleChargedCut05GeV(G4double theta)            {if (fBoolChargedCut05GeV == 0)   { fAngleChargedCut05GeV  += theta; fBoolChargedCut05GeV = 1; } }; // If angle < 30 mrad
    void AddAngleChargedCut1GeV(G4double theta)             {if (fBoolChargedCut1GeV == 0)   { fAngleChargedCut1GeV  += theta; fBoolChargedCut1GeV = 1; } }; // If angle < 30 mrad
    
   
    void AddAngleIoni(G4double theta)             {if (fBooleIoni == 0)   { fAngleIoni  += theta; fBooleIoni = 1; } }; // If angle < 30 mrad
    void AddAnglePP(G4double theta)             {
       if ((fBoolePP1 == 0) && (fNbTracks == 1))   { fAnglePP1  += theta; fBoolePP1 = 1; }
       if ((fBoolePP2 == 0) && (fNbTracks == 2))   { fAnglePP2  += theta; fBoolePP2 = 1; }
    }; // If angle < 30 mrad 
    
    void AddAngleIoniCut05GeV(G4double theta)             {if (fBooleIoniCut05GeV == 0)   { fAngleIoniCut05GeV  += theta; fBooleIoniCut05GeV = 1; } }; // If angle < 30 mrad
    void AddAnglePPCut05GeV(G4double theta)             {
       if ((fBoolePP1Cut05GeV == 0) && (fNbTracksCut05GeV == 1))   { fAnglePP1Cut05GeV  += theta; fBoolePP1Cut05GeV = 1; }
       if ((fBoolePP2Cut05GeV == 0) && (fNbTracksCut05GeV == 2))   { fAnglePP2Cut05GeV  += theta; fBoolePP2Cut05GeV = 1; }
    }; // If angle < 30 mrad 
    
    void AddAngleIoniCut1GeV(G4double theta)             {if (fBooleIoniCut1GeV == 0)   { fAngleIoniCut1GeV  += theta; fBooleIoniCut1GeV = 1; } }; // If angle < 30 mrad
    void AddAnglePPCut1GeV(G4double theta)             {
       if ((fBoolePP1Cut1GeV == 0) && (fNbTracksCut1GeV == 1))   { fAnglePP1Cut1GeV  += theta; fBoolePP1Cut1GeV = 1; }
       if ((fBoolePP2Cut1GeV == 0) && (fNbTracksCut1GeV == 2))   { fAnglePP2Cut1GeV  += theta; fBoolePP2Cut1GeV = 1; }
    }; // If angle < 30 mrad 
       
           
    //for separate processes
    void AddNbTracks()                                      {fNbTracks++;}  // If angles < 30 mrad
    void AddNbTracksCut05GeV()                              {fNbTracksCut05GeV++;}  // If angles < 30 mrad
    void AddNbTracksCut1GeV()                               {fNbTracksCut1GeV++;}  // If angles < 30 mrad
    

    /*void AddMomentumDirectionPrimExT(G4double dirx, G4double diry, G4double dirz)  {fMomDirPx += dirx; fMomDirPy += diry; fMomDirPz += dirz;};
    void AddMomentumDirectionSecExT(G4double dirx, G4double diry, G4double dirz)  {fMomDirSx += dirx; fMomDirSy += diry; fMomDirSz += dirz;};*/
       
  private:
    DetectorConstruction* fDetectorconstruction;
    RunAction*    fRunAction;
    
    G4double      fEnergyExitTargetPrimary;
    G4double      fEnergyExitTargetSecondaryPositron;
    G4double      fEnergyExitTargetSecondaryElectron;
    G4double      fEnergyExitTargetSecondaryMuonPlus;
    G4double      fEnergyExitTargetSecondaryMuonMinus;
    G4double      fEnergyExitTargetCharged;
    G4double      fEnergyExitTargetChargedCut05GeV;
    G4double      fEnergyExitTargetChargedCut1GeV;
    
    G4double      fEnergyExitTargetIoni;
    G4double      fEnergyExitTargetPP1;
    G4double      fEnergyExitTargetPP2;
    
    G4double      fEnergyExitTargetIoniCut05GeV;
    G4double      fEnergyExitTargetPP1Cut05GeV;
    G4double      fEnergyExitTargetPP2Cut05GeV;
    
    G4double      fEnergyExitTargetIoniCut1GeV;
    G4double      fEnergyExitTargetPP1Cut1GeV;
    G4double      fEnergyExitTargetPP2Cut1GeV;
    
    G4double      fAnglePrimary;
    G4double      fAngleSecondaryPositron;
    G4double      fAngleSecondaryElectron;
    G4double      fAngleSecondaryMuonPlus;
    G4double      fAngleSecondaryMuonMinus;
    G4double      fAngleCharged;
    G4double      fAngleChargedCut05GeV;
    G4double      fAngleChargedCut1GeV;
    
    G4double      fAngleIoni;
    G4double      fAnglePP1;
    G4double      fAnglePP2;
    
    G4double      fAngleIoniCut05GeV;
    G4double      fAnglePP1Cut05GeV;
    G4double      fAnglePP2Cut05GeV;
    
    G4double      fAngleIoniCut1GeV;
    G4double      fAnglePP1Cut1GeV;
    G4double      fAnglePP2Cut1GeV;

    G4int         fBoolPrim, fBoolSecPositron, fBoolSecElectron, fBoolSecMuonPlus, fBoolSecMuonMinus, fBoolCharged, fBoolChargedCut05GeV, fBoolChargedCut1GeV, fBooleIoni, fBoolePP1, fBoolePP2, fBooleIoniCut05GeV, fBoolePP1Cut05GeV, fBoolePP2Cut05GeV, fBooleIoniCut1GeV, fBoolePP1Cut1GeV, fBoolePP2Cut1GeV; //a boolean which becomes 1 when a primary or secondary, respectively, is found with angle < 30 mrad
    G4int         fNbTracks, fNbTracksCut05GeV, fNbTracksCut1GeV;

    /*G4double      fMomDirPx, fMomDirPy, fMomDirPz;
    G4double      fMomDirSx, fMomDirSy, fMomDirSz;*/
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
