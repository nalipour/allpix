/**
 *  Author:
 *    
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixMCTruthDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

AllPixMCTruthDigitizer::AllPixMCTruthDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
  : AllPixDigitizerInterface (modName) {

  // Registration of digits collection name
  collectionName.push_back(digitColName);
  m_hitsColName.push_back(hitsColName);

  // threshold
  m_digitIn.thl = 0.;

}

AllPixMCTruthDigitizer::~AllPixMCTruthDigitizer(){

}

void AllPixMCTruthDigitizer::Digitize(){

  // create the digits collection
  m_digitsCollection = new AllPixMCTruthDigitsCollection("AllPixMCTruthDigitizer", collectionName[0] );

  // get the digiManager
  G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

  // BoxSD_0_HitsCollection
  G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

  AllPixTrackerHitsCollection * hitsCollection = 0;
  hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

  // temporary data structure
  map<pair<G4int, G4int>, G4double > pixelsContent;
  pair<G4int, G4int> tempPixel;

  G4int nEntries = hitsCollection->entries();

  // Example of detector description handle
  // provided by the interface
  AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
  gD->GetNPixelsX();
  double thickness=gD->GetSensorZ(); //thickness[mm] 


  G4double AvgPosX=0.0;
  G4double AvgPosY=0.0;
  G4double AvgPosZ=0.0;

  G4double LocalAvgPosX=0.0;
  G4double LocalAvgPosY=0.0;
  G4double LocalAvgPosZ=0.0;

  G4double x0, y0, z0, xp, yp, zp=0;

  G4double MC_deposited_energy=0.0;
  int pixelX=-10;
  int pixelY=-10;

  for(G4int itr  = 0 ; itr < nEntries ; itr++) {

    tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
    tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
    pixelsContent[tempPixel] += (*hitsCollection)[itr]->GetEdep();

    MC_deposited_energy+=(*hitsCollection)[itr]->GetEdep();

    G4double xpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().x();
    G4double ypos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().y();
    G4double zpos=(*hitsCollection)[itr]->GetPosWithRespectToPixel().z()+thickness/2.0; // [mm]; zpos=thickness corresponds to the sensor side and zpos=0 corresponds to the pixel side

    G4double global_xpos=(*hitsCollection)[itr]->GetGlobalPosition().x();
    G4double global_ypos=(*hitsCollection)[itr]->GetGlobalPosition().y();
    G4double global_zpos=(*hitsCollection)[itr]->GetGlobalPosition().z();


      G4double local_xpos=(*hitsCollection)[itr]->GetLocalPosition().x();
      G4double local_ypos=(*hitsCollection)[itr]->GetLocalPosition().y();
      G4double local_zpos=(*hitsCollection)[itr]->GetLocalPosition().z();

      pixelX=(*hitsCollection)[itr]->GetPixelNbX();
      pixelY=(*hitsCollection)[itr]->GetPixelNbY();

    if (itr==0)
      {
	AvgPosX=global_xpos;
	AvgPosY=global_ypos;
	AvgPosZ=global_zpos;

	LocalAvgPosX=local_xpos;
	LocalAvgPosY=local_ypos;
	LocalAvgPosZ=local_zpos;
      }
  }

  map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

  // NOTE that there is a nice interface which provides useful info for hits.
  // For instance, the following member gives you the position of a hit with respect
  //  to the center of the pixel.
  // G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
  // See class AllPixTrackerHit !

  // Also, you have full access to the detector geometry from this scope
  // AllPixGeoDsc * GetDetectorGeoDscPtr()
  // provides you with a pointer to the geometry description.
  // See class AllPixGeoDsc !

  // for( ; pCItr != pixelsContent.end() ; pCItr++)
  //   {

  //     if((*pCItr).second > 0) // over threshold !
  // 	{
  // 	  // Create one digit per pixel, I need to look at all the pixels first
  // 	  AllPixMCTruthDigit * digit = new AllPixMCTruthDigit;
  // 	  digit->SetPixelIDX((*pCItr).first.first);
  // 	  digit->SetPixelIDY((*pCItr).first.second);
  // 	  digit->SetPixelCounts((*pCItr).second/keV);
  // 	  G4cout << "Energy=" << (*pCItr).second/keV << G4endl;

  // 	  G4cout << "X: " << (*pCItr).first.first << ", Y: " << (*pCItr).first.second << G4endl;
  // 	  G4cout << "Global: x=" << AvgPosX << ", y=" << AvgPosY << ", z=" << AvgPosZ << G4endl;
  // 	  G4cout << "Local: x=" << LocalAvgPosX << ", y=" << LocalAvgPosY << ", z=" << LocalAvgPosZ << G4endl;


  // 	  // ====== TO be corrected later==================== //
  // 	  digit->Set_posX_WithRespectoToPixel(AvgPosX);// /nEntries);
  // 	  digit->Set_posY_WithRespectoToPixel(AvgPosY); // /nEntries);
  // 	  digit->Set_posZ_WithRespectoToPixel(AvgPosZ); // /nEntries);
  // 	  //===================================================//
  // 	  // G4cout << "dEdX : " << (*pCItr).second/keV << endl;

  // 	  m_digitsCollection->insert(digit);
  // 	}
  //   }

  if(MC_deposited_energy>0)
    {
      AllPixMCTruthDigit * digit = new AllPixMCTruthDigit;
      digit->SetPixelEnergyDep(MC_deposited_energy/keV);
      digit->SetPixelIDX(pixelX); // An ugly hack to have the same coordinate system as the test-beam
      digit->SetPixelIDY(pixelY);	 
      digit->IncreasePixelCounts(); // Counting mode
      m_digitsCollection->insert(digit);
	
    }

  G4int dc_entries = m_digitsCollection->entries();
  if(dc_entries > 0){
    G4cout << "--------> Digits Collection : " << collectionName[0]
	   << "(" << m_hitsColName[0] << ")"
	   << " contains " << dc_entries
	   << " digits" << G4endl;
  }

  StoreDigiCollection(m_digitsCollection);

}
