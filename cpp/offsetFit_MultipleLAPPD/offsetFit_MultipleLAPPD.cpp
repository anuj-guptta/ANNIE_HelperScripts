/*
==================================================
Script Name     : offsetFit_MultipleLAPPD.cpp
Author          : Yue Feng
Created On      : N/A
Updated By      : Anuj Gupta
Last Updated    : 2026-09-19

Purpose         : The script performs LAPPD timing offset fitting and timing corrections.
                  Detailed procedure is described in README.md

Inputs          : It takes 5 arguments as input.
                  filename: LAPPDTree.root, fitTargetTriggerWord, triggerGrouped, intervalInSeconds, processPfNumber

Outputs         : Produces offsetFitResult.root (used during Event Building in LAPPDLoadStore tool)
                  Other outputs: offsetFit.txt, outputEvents.txt

Usage           : root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 14, 1, 10, 0)' // for beam runs
                  root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 47, 0, 10, 0)' // for laser runs
                  Please refer to README.md for additional details.

Dependencies    : ROOT
==================================================
*/

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"

vector<vector<ULong64_t>> fitInThisReset(
    const std::vector<ULong64_t> &LAPPDDataTimeStampUL,
    const std::vector<ULong64_t> &LAPPDDataBeamgateUL,
    const std::vector<ULong64_t> &LAPPD_PPS,
    const int fitTargetTriggerWord,
    const std::vector<ULong64_t> &CTCTrigger,
    const std::vector<ULong64_t> &CTCPPS,
    const ULong64_t PPSDeltaT, // here, pre-scale is in units ps
    const int partFileNumber,
    const int LAPPD_ID,
    const int ACDCNumber)
{
    std::cout << "***************************************" << endl;
    std::cout << "Fitting part file number: " << partFileNumber << ", ACDC number: " << ACDCNumber << std::endl;
    std::cout << "Fitting in this reset with:" << std::endl;
    std::cout << "LAPPDDataTimeStampUL size: " << LAPPDDataTimeStampUL.size() << std::endl;
    std::cout << "LAPPDDataBeamgateUL size: " << LAPPDDataBeamgateUL.size() << std::endl;
    std::cout << "LAPPD_PPS size: " << LAPPD_PPS.size() << std::endl;
    std::cout << "CTCTrigger size: " << CTCTrigger.size() << std::endl;
    std::cout << "CTCPPS size: " << CTCPPS.size() << std::endl;
    std::cout << "PPSDeltaT (us): " << PPSDeltaT / 1E6 << std::endl;
    
    // For this input PPS, fit an offset
    // Return the offset and other information in order
    // Precedure:
    	// 1. Check the drift
    	// 2. Shift the timestamp and beamgate based on drift
    	// 3. Fit the offset

    std::vector<double> PPSInterval_ACDC;
    for (size_t i = 1; i < LAPPD_PPS.size(); i++)
    {
	// This condition is already checked before calling this function.
	// I am adding this condition here again for debugging purpose.
	if (LAPPD_PPS[i] < LAPPD_PPS[i - 1]) {
		std::cout << "WARNING: unexpected PPS ordering violation in "
             		<< "part file " << partFileNumber
             		<< ", ACDC " << ACDCNumber
             		<< ", index " << i << std::endl;
		continue;
	}

        double diff = static_cast<double>(LAPPD_PPS[i] - LAPPD_PPS[i - 1]);
	std::cout << fixed << "PPSDiff (us) = " << diff / 1E6 << std::endl;
        
	// only save the time interval that is close to the deltaT, to avoid the clock tick level miss and the whole PPS pulse was missed.
        if (diff > 0.9 * PPSDeltaT && diff < 1.1 * PPSDeltaT)
        {
            PPSInterval_ACDC.push_back((PPSDeltaT - diff) / 1E6); // in microseconds
            std::cout << "push back interval (us) = " << (PPSDeltaT - diff) / 1E6 
		    << ", 0.9*PPSdeltaT  = " << 0.9 * PPSDeltaT / 1E6 
		    << ", 1.1*PPSdeltaT = " << 1.1 * PPSDeltaT / 1E6 << std::endl;
        }
        else
        {
		std::cout << "In units of us, 0.9*PPSdeltaT = " << 0.9 * PPSDeltaT / 1E6 
		    << ", 1.1*PPSdeltaT = " << 1.1 * PPSDeltaT / 1E6
		    << ", PPSDeltaT - diff = " << (PPSDeltaT - diff) / 1E6 << std::endl;
        }
    }
   
    // Create a mathematically UNIQUE name for the histogram
    TString hName = TString::Format("h_PF%d_LAPPD%d_ACDC%d", partFileNumber, LAPPD_ID, ACDCNumber);
    TH1D *h = new TH1D(hName, hName, 1000, 0, 1E3);

    for (size_t i = 0; i < PPSInterval_ACDC.size(); i++)
    {
        // only fill the drift histogram if there is a drift > 2 microseconds
	std::cout << "PPSInterval_ACDC[" << i << "] = " << PPSInterval_ACDC[i] << std::endl;
        if (PPSInterval_ACDC[i] > 2) 
        {
            h->Fill(PPSInterval_ACDC[i]); // fill in microseconds
	    std::cout << "Fill histogram: PPSInterval_ACDC[" << i << "]: " << PPSInterval_ACDC[i] << std::endl;
        }
    }

    // Create a mathematically UNIQUE name for the fit function
    TString fName = TString::Format("gausf_PF%d_LAPPD%d_ACDC%d", partFileNumber, LAPPD_ID, ACDCNumber);
    TF1 *gausf = new TF1(fName, "gaus", 0, 1E3);

    // Perform the fit
    h->Fit(gausf, "Q"); // Q for quiet mode
    
    ULong64_t drift = static_cast<ULong64_t>(gausf->GetParameter(1) * 1E6); // drift in ps
    ULong64_t trueInterval = PPSDeltaT - drift;
    std::cout << "Gaussian Drift in ps is " << drift << std::endl;
    std::cout << "True PPS interval in ps is " << trueInterval << std::endl;
    
    // Delete objects to avoid memory leak  
    delete gausf;
    delete h;

    std::map<double, std::vector<vector<int>>> DerivationMap;
    std::cout << "LAPPD_PPS.size() " << LAPPD_PPS.size() << ", CTCPPS.size() " << CTCPPS.size() << std::endl;
    
    for (int i = 0; i < LAPPD_PPS.size(); i++)
    {
        if (i > 5)
        {
            if (i % static_cast<int>(LAPPD_PPS.size() / 5) == 0)
                cout << "Fitting PPS " << i << " of " << LAPPD_PPS.size() << endl;
        }

        ULong64_t LAPPD_PPS_ns = LAPPD_PPS.at(i) / 1000; // ns units
        ULong64_t LAPPD_PPS_truncated_ps = LAPPD_PPS.at(i) % 1000;

	// CTCPPS are in ns units
        for (int j = 0; j < CTCPPS.size(); j++)
        {
            vector<double> diffSum;
            ULong64_t offsetNow_ns = 0;

            if (drift == 0)
            {
                offsetNow_ns = CTCPPS.at(j) - LAPPD_PPS_ns;
            }
            else
            {
                double LAPPD_PPS_ns_dd = static_cast<double>(LAPPD_PPS_ns);
                double driftScaling = LAPPD_PPS_ns_dd / (trueInterval * 1.0 / 1000);
                ULong64_t totalDriftedClock = static_cast<ULong64_t>(static_cast<double>(drift) * driftScaling / 1000);
                offsetNow_ns = CTCPPS.at(j) - (LAPPD_PPS_ns + totalDriftedClock);
                
		if (i < 3 && j < 3) {
		    std::cout << "LAPPD_PPS_ns = " << LAPPD_PPS_ns
     		    	<< ", trueInterval(ns) = " << trueInterval/1000
     			<< ", elapsed_interval = " << driftScaling
     			<< ", drift_per_interval(ns) = " << drift/1000
     			<< ", totalDriftedClock(ns) = " << totalDriftedClock << std::endl;
		}
	    }

            if (i < 3 && j < 3) {
                std::cout << "In ns: CTCPPS.at(" << j << "): " << CTCPPS.at(j) 
			<< ", LAPPD_PPS.at(" << i << "): " << LAPPD_PPS.at(i) 
			<< ", use LAPPD_PPS_ns: " << LAPPD_PPS_ns 
			<< ", offsetNow_ns: " << offsetNow_ns << std::endl;
	    }

            // now, we have offset in ns.
            // We do fit in ns level, then save the precise information in ps level.

   	    // count the number that how many data event timestamp doesn't matched to a target trigger within an interval.
            int orphanCount = 0;
            vector<int> notOrphanIndex;
            vector<int> ctcPairedIndex;
            vector<int> ctcOrphanPairedIndex;

            for (int lappdb = 0; lappdb < LAPPDDataBeamgateUL.size(); lappdb++)
            {
                ULong64_t TS_ns = LAPPDDataTimeStampUL.at(lappdb) / 1000;
                ULong64_t TS_truncated_ps = LAPPDDataTimeStampUL.at(lappdb) % 1000;
                
		// if fit for the undelayed beam trigger (14), use Beamgate
                if (fitTargetTriggerWord == 14)
                {
                    TS_ns = LAPPDDataBeamgateUL.at(lappdb) / 1000;
                    TS_truncated_ps = LAPPDDataBeamgateUL.at(lappdb) % 1000;
                }

                ULong64_t driftCorrectionForTS = 0;
                if (drift != 0)
                {
                    double TS_ns_dd = static_cast<double>(TS_ns);
                    double driftScaling = TS_ns_dd / (trueInterval * 1.0 / 1000);
                    driftCorrectionForTS = static_cast<ULong64_t>(static_cast<double>(drift) * driftScaling / 1000);
                    
		    if (lappdb<3 && i<3 && j<3) {
			    std::cout << "i = " << i
     				<< ", j = " << j
     				<< ", lappdb = " << lappdb
     				<< ", elapsed_interval = " << driftScaling
     				<< ", drift_per_interval(ns) = " << drift/1000
     				<< ", driftCorrectionForTS(ns) = " << driftCorrectionForTS << std::endl;
		    }
		}

                // this is the TS we use for matching
                ULong64_t DriftCorrectedTS_ns = TS_ns + offsetNow_ns + driftCorrectionForTS;
		if (lappdb < 3 && i < 3 && j < 3)
		{
			std::cout << "LAPPD Timestamp used for matching: " << DriftCorrectedTS_ns << std::endl;
		}

                // initialise other matching variables
                Long64_t minMatchDiff = std::numeric_limits<Long64_t>::max();
                int minPairIndex = 0;
                string reason = "none";
                Long64_t firstDiff = 0;
                bool useValue = true;

                // find the best match of target trigger to this TS, in ns level.
                for (int ctcb = 0; ctcb < CTCTrigger.size(); ctcb++)
                {
                    ULong64_t CTCTrigger_ns = CTCTrigger.at(ctcb);
                    Long64_t diff = static_cast<Long64_t>(CTCTrigger_ns) - static_cast<Long64_t>(DriftCorrectedTS_ns);
                    
		    if (diff < 0)
                        diff = -diff;

                    if (diff < minMatchDiff)
                    {
                        minMatchDiff = diff;
                        minPairIndex = ctcb;
                    }

                    Long64_t LastBound = diff - minMatchDiff;
                    if (LastBound < 0)
                        LastBound = -LastBound;
                    Long64_t FirstBound = diff - firstDiff;
                    if (FirstBound < 0)
                        FirstBound = -FirstBound;

                    // save the dt for the first matching to determin the matching position in range.
                    if (ctcb == 0)
                        firstDiff = diff;

                    if (ctcb == CTCTrigger.size() - 1 && LastBound / 1E9 < 0.01)
                    {
                        useValue = false;
                        reason = "When matching the last TS, this diff is too close to (or even is) the minMatchDiff at index " 
				+ std::to_string(minPairIndex) + " in " + std::to_string(CTCTrigger.size()) 
				+ ". All LAPPD TS is out of the end of CTC trigger range";
                    }
                    if (ctcb == CTCTrigger.size() - 1 && FirstBound / 1E9 < 0.01)
                    {
                        useValue = false;
                        reason = "When matching the last TS, this diff is too close to (or even is) the firstDiff at index " 
				+ std::to_string(minPairIndex) + " in " + std::to_string(CTCTrigger.size()) 
				+ ". All LAPPD TS is out of the beginning of CTC trigger range";
                    }
                }

                if (useValue)
                {
                    diffSum.push_back(minMatchDiff);
                    // std::cout << "LAPPD PPS " << i << " = " << LAPPD_PPS.at(i) 
		    	// << " ps, CTC PPS " << j << " = " << CTCPPS.at(j)  << " ns , LAPPDDataBeamgateUL[" << lappdb << "] = " << LAPPDDataBeamgateUL[lappdb]
			// <<" ps, CTCTrigger[" << minPairIndex << "] = " << CTCTrigger[minPairIndex] << " ns, minMatchDiff = " << minMatchDiff << " ns" << std::endl;
                    
		    double minAllowedDiff = 0;
                    double maxAllowedDiff = 100E3;
                    if (fitTargetTriggerWord == 14)
                    {
                        minAllowedDiff = 322E3;
                        maxAllowedDiff = 326E3;
                    }
                    if (minMatchDiff > maxAllowedDiff || minMatchDiff < minAllowedDiff)
                    {
                        // TODO: adjust the limit for laser.
                        orphanCount += 1;
                        ctcOrphanPairedIndex.push_back(minPairIndex);
                    }
                    else
                    {
			// std::cout << minMatchDiff << std::endl;
                        // notOrphanIndex.push_back(lappdb);
                        notOrphanIndex.push_back(diffSum.size() - 1);
			ctcPairedIndex.push_back(minPairIndex);
                    }
                }
                else
                {
                    orphanCount += 1;
                }
            }

	    // std::cout << "i = " << i << " | j = " << j << " | diffSum = " << diffSum.size() << " | orphanCount = " << orphanCount << std::endl;

            double mean_dev = 0;
            for (int k = 0; k < diffSum.size(); k++)
            {
                mean_dev += diffSum[k];
		// std::cout << "DiffSum[" << k << "] = " << diffSum[k] << std::endl;
            }
            if (diffSum.size() > 0)
                mean_dev = mean_dev / diffSum.size();
        
            double mean_dev_noOrphan = 0;
            int pairedCount = 0;
            for (int k = 0; k < notOrphanIndex.size(); k++)
            {
		int idx = notOrphanIndex[k];
		if (idx < 0 || idx >= diffSum.size())
    		{
       	 		std::cout << "SKIP BAD IDX: " << idx << std::endl;
        		continue;
    		}
		
		if (fitTargetTriggerWord == 14)
    		{
                        if (diffSum[idx] > 322E3 && diffSum[idx] < 326E3)
        		{
            		
				mean_dev_noOrphan += diffSum[idx];
            			pairedCount += 1;
        		}
    		}
		else {
			mean_dev_noOrphan += diffSum[idx];
			pairedCount += 1;
		}
            }
            
            if (pairedCount != 0)
            {
                mean_dev_noOrphan = mean_dev_noOrphan / pairedCount;
            }
            else
            {
                mean_dev_noOrphan = -1; // if all timestamps are out of the range, set it to -1
            }

            const int maxAttempts = 1000;
            int attemptCount = 0;

            // The mean_dev can't be larger than 1 s because the beam spill is 15Hz.
            // So, a good match will give a large number of events in desired range on integer level, minus the mean _dev.
            // By selecting the largest quality Number, for matchs with the same matched beamgate, smaller mean_dev is better.
            double qualityNumber = pairedCount*1e9 - mean_dev;

            bool debug_print = false;
            if (debug_print)
            {            
		// print all diffSum
		std::cout << "diffSum: ";
                for (const auto &diff : diffSum)
                {
			std::cout << diff << ", ";
                }
                std::cout << std::endl;
            }

            if (mean_dev > 0)
            {
                int increment_dev = 0;
                while (attemptCount < maxAttempts)
                {
                    // TODO
                    // auto iter = DerivationMap.find(mean_dev);
                    auto iter = DerivationMap.find(qualityNumber);
		    if (iter == DerivationMap.end() || iter->second.empty())
                    {
                        vector<int> Info = {i, j, orphanCount, static_cast<int>(mean_dev_noOrphan * 1000), increment_dev};
                        /*
                        DerivationMap[mean_dev].push_back(Info);
                        DerivationMap[mean_dev].push_back(notOrphanIndex);
                        DerivationMap[mean_dev].push_back(ctcPairedIndex);
                        DerivationMap[mean_dev].push_back(ctcOrphanPairedIndex);
                        */

                        DerivationMap[qualityNumber].push_back(Info);
                        DerivationMap[qualityNumber].push_back(notOrphanIndex);
                        DerivationMap[qualityNumber].push_back(ctcPairedIndex);
                        DerivationMap[qualityNumber].push_back(ctcOrphanPairedIndex);

                        break;
                    }
                    else
                    {
                        increment_dev += 1;
                        // mean_dev += 0.001; // if the mean_dev is already in the map, increase it by 1ps
                        qualityNumber += 1e-6;
                        attemptCount += 1;
                    }
                }
            }
        }
    }

    // finish matching, found the maximum qualityNumber in the map, extract the matching information
    double min_mean_dev = std::numeric_limits<double>::max();
    double max_qualityNumber = 0;
    int final_i = 0;
    int final_j = 0;
    int gotOrphanCount = 0;
    double gotMin_mean_dev_noOrphan = 0;
    double increament_times = 0;
    vector<int> final_notOrphanIndex;
    vector<int> final_ctcPairedIndex;
    vector<int> final_ctcOrphanIndex;
    
    for (const auto &minIter : DerivationMap)
    {
        //if (minIter.first > 10 && minIter.first < min_mean_dev)
        if (minIter.first > -1 && minIter.first > max_qualityNumber)
        {
            //min_mean_dev = minIter.first;
            max_qualityNumber = minIter.first;
            double qnum_ns = minIter.first / 1e9;
            min_mean_dev = (std::ceil(qnum_ns) - qnum_ns)*1e9;
            
	    final_i = minIter.second[0][0];
            final_j = minIter.second[0][1];
            gotOrphanCount = minIter.second[0][2];
            gotMin_mean_dev_noOrphan = static_cast<double>(minIter.second[0][3]) / 1000;
            increament_times = minIter.second[0][4];
            
	    final_notOrphanIndex = minIter.second[1];
            final_ctcPairedIndex = minIter.second[2];
            final_ctcOrphanIndex = minIter.second[3];
        }
    }

    ULong64_t final_offset_ns = 0;
    ULong64_t final_offset_ps_negative = 0;
    if (drift == 0)
    {
	std::cout << "Drift is 0 above ns level." << " Using CTC PPS at index " << final_j 
		<< " = " << CTCPPS.at(final_j) << " and LAPPD PPS at index " << final_i << " = " << LAPPD_PPS.at(final_i)/1000 << std::endl;
        
	final_offset_ns = CTCPPS.at(final_j) - (LAPPD_PPS.at(final_i) / 1000);
        final_offset_ps_negative = LAPPD_PPS.at(final_i) % 1000;
    }
    else
    {
        ULong64_t LAPPD_PPS_ns = LAPPD_PPS.at(final_i) / 1000;
        ULong64_t LAPPD_PPS_truncated_ps = LAPPD_PPS.at(final_i) % 1000;
        double driftScaling = static_cast<double>(LAPPD_PPS_ns) / (trueInterval * 1.0 / 1000); // this is the same drift scaling as in the matching loop
        ULong64_t totalDriftedClock = static_cast<ULong64_t>(static_cast<double>(drift) * driftScaling / 1000);
        std::cout << fixed << "Found CTCPPS as " << CTCPPS.at(final_j) << ", LAPPD_PPS_ns = " << LAPPD_PPS_ns;

        final_offset_ns = CTCPPS.at(final_j) - (LAPPD_PPS_ns + totalDriftedClock);
        final_offset_ps_negative = LAPPD_PPS_truncated_ps; // this is useless if there is a drift though.
	
	std::cout << " final_offset_ns = " << final_offset_ns << ", driftScaling = "
		<< driftScaling << ", totalDriftedClock_ns = " << totalDriftedClock << std::endl;
    }

    /*
    cout << "******* Fit Finished *******" << endl;
    cout << "*** Final offset in is " << final_offset_ns << " ns minus " << final_offset_ps_negative << " ps" << endl;
    cout << "*** Final orphan count is " << gotOrphanCount << endl;
    cout << "*** Final mean_dev_noOrphan is " << gotMin_mean_dev_noOrphan << " ps" << endl;
    cout << "*** Final increament times in this result is " << increament_times << endl;
    cout << "*** Final mean deviation is " << min_mean_dev << " ps" << endl;
    cout << "*** Final PPS index is " << final_i << ", in total of " << LAPPD_PPS.size() << endl;
    cout << "*** Final CTC PPS index is " << final_j << ", in total of " << CTCPPS.size() << endl;
    cout << "***************************" << endl;
    cout << "******* Saving *************" << endl;
    */

    cout << "\033[1;34m******* Fit Finished *******\033[0m" << endl;
    cout << "\033[1;34m*** Part file number is \033[1;31m" << partFileNumber << "\033[1;34m, ACDC number is \033[1;31m" << ACDCNumber << "\033[0m" << endl;
    cout << "\033[1;34m*** Final offset in is \033[1;31m" << final_offset_ns << "\033[1;34m ns minus \033[1;31m" << final_offset_ps_negative << "\033[1;34m ps\033[0m" << endl;
    cout << "\033[1;34m*** Final orphan count is \033[1;31m" << gotOrphanCount << "\033[0m" << endl;
    cout << "\033[1;34m*** Final mean_dev_noOrphan is \033[1;31m" << gotMin_mean_dev_noOrphan << "\033[1;34m ns\033[0m" << endl;
    cout << "\033[1;34m*** Final increament times in this result is \033[1;31m" << increament_times << "\033[0m" << endl;
    cout << "\033[1;34m*** Final quality number is \033[1;31m" << max_qualityNumber << "\033[0m" << endl;
    cout << "\033[1;34m*** Final mean deviation is \033[1;31m" << min_mean_dev << "\033[1;34m ns\033[0m" << endl;
    cout << "\033[1;34m*** Final PPS index is \033[1;31m" << final_i << "\033[1;34m, in total of \033[1;31m" << LAPPD_PPS.size() << "\033[0m" << endl;
    cout << "\033[1;34m*** Final CTC PPS index is \033[1;31m" << final_j << "\033[1;34m, in total of \033[1;31m" << CTCPPS.size() << "\033[0m" << endl;
    
    std::cout << "ACDC " << ACDCNumber << " matched CTC triggers: ";
    for (auto x : final_ctcPairedIndex)
        std::cout << x << " ";

    std::cout << std::endl;
    std::cout << "Final_i = " << final_i
          << "   LAPPD_PPS = " << LAPPD_PPS.at(final_i)
          << std::endl;

    std::cout << "Final_j = " << final_j
          << "   CTC_PPS = " << CTCPPS.at(final_j)
          << std::endl;

    cout << "\033[1;34m***************************\033[0m" << endl;
    cout << "\033[1;34m******* Saving *************\033[0m" << endl;

    // now, based on this offset, calculate the event time for each event, and save to result.
    vector<ULong64_t> TimeStampRaw;
    vector<ULong64_t> BeamGateRaw;
    vector<ULong64_t> TimeStamp_ns;
    vector<ULong64_t> BeamGate_ns;
    vector<ULong64_t> TimeStamp_ps;
    vector<ULong64_t> BeamGate_ps;
    vector<ULong64_t> EventIndex;
    vector<ULong64_t> EventDeviation_ns;
    vector<ULong64_t> CTCTriggerIndex;
    vector<ULong64_t> CTCTriggerTimeStamp_ns;
    vector<ULong64_t> BeamGate_correction_tick;
    vector<ULong64_t> TimeStamp_correction_tick;
    vector<long long> PPS_tick_correction; // if timestamp falls in between of PPS i and i+1, corrected timestamp = timestamp + PPS_tick_correction[i]
    vector<long long> LAPPD_PPS_missing_ticks;
    vector<ULong64_t> LAPPD_PPS_interval_ticks;

    vector<ULong64_t> BG_PPSBefore;
    vector<ULong64_t> BG_PPSAfter;
    vector<ULong64_t> BG_PPSDiff;
    vector<ULong64_t> BG_PPSMiss;
    vector<ULong64_t> TS_PPSBefore;
    vector<ULong64_t> TS_PPSAfter;
    vector<ULong64_t> TS_PPSDiff;
    vector<ULong64_t> TS_PPSMiss;

    vector<ULong64_t> TS_driftCorrection_ns;
    vector<ULong64_t> BG_driftCorrection_ns;

    // calculate the missing ticks for each LAPPD PPS
    PPS_tick_correction.push_back(0);
    LAPPD_PPS_missing_ticks.push_back(0);
    LAPPD_PPS_interval_ticks.push_back(0);
    ULong64_t intervalTicks = 320000000 * (PPSDeltaT / 1E12);
    const ULong64_t LAPPD_PPS_SCALE = 3125;

    for (int i = 1; i < LAPPD_PPS.size(); i++)
    {
        ULong64_t thisPPS = LAPPD_PPS.at(i) / LAPPD_PPS_SCALE;
        ULong64_t prevPPS = LAPPD_PPS.at(i - 1) / LAPPD_PPS_SCALE;
        
	// std::cout << "LAPPD_PPS.at(" << i << ") = " << LAPPD_PPS.at(i) << ", LAPPD_PPS.at(" << i-1 << ") = " << LAPPD_PPS.at(i-1) << std::endl;
	long long thisInterval = static_cast<long long>(thisPPS) 
		- static_cast<long long>(prevPPS);
        long long thisMissingTicks = static_cast<long long>(intervalTicks) - thisInterval;
        LAPPD_PPS_interval_ticks.push_back(thisInterval);

	// if the difference between this PPS and previous PPS is > intervalTicks-20 and < intervalTicks+5,
        // the missed value is thisInterval - intervalTicks
        // push this value plus the sum of previous missing ticks to the vector
        // else just push the sum of previous missing ticks
	
        long long sumOfPreviousMissingTicks = 0;
        for (int j = 0; j < LAPPD_PPS_missing_ticks.size(); j++)
        {
            sumOfPreviousMissingTicks += LAPPD_PPS_missing_ticks.at(j);
        }
        
        if (thisMissingTicks > -20 && thisMissingTicks < 20)
        {
            LAPPD_PPS_missing_ticks.push_back(thisMissingTicks);
            PPS_tick_correction.push_back(thisMissingTicks + sumOfPreviousMissingTicks);
        }
        else
        {
            // some time one pps might be wired, but the combination of two is ok.
            bool combined = false;
            
	    // if there are missing PPS, interval is like 22399999990 % 3200000000 = 3199999990
            long long pInterval = thisInterval % static_cast<long long>(intervalTicks);
            long long missingPTicks = 0;
           
	    // std::cout << "pInterval = " << pInterval << std::endl;
	   
	    // Case 1: There is a missing heartbeat and is an exact multiple of IntervalTicks
	    if (pInterval == 0)
	    {
	        LAPPD_PPS_missing_ticks.push_back(0);
	        PPS_tick_correction.push_back(sumOfPreviousMissingTicks);
		std::cout << "Pushing PPS correction " << i << ", this PPS interval tick is " << thisInterval 
			<< ", missing ticks: " << thisMissingTicks << ", push missing " << LAPPD_PPS_missing_ticks.back() 
			<< ", push correction " << PPS_tick_correction.back() << std::endl;
		continue;
	    }

	    // Case 2: There is a missing heartbeat but the next signal has a small drift
	    if (static_cast<long long>(intervalTicks) - pInterval < 30)
                missingPTicks = static_cast<long long>(intervalTicks) - pInterval;
            else if (pInterval < 30)
                missingPTicks = -pInterval;

	    // std::cout << "missingPTicks = " << missingPTicks << std::endl;
            //////
            if (missingPTicks == 1 || missingPTicks == -1)
            {
		std::cout << "Found missing tick is " << missingPTicks << ", continue." << std::endl;
                LAPPD_PPS_missing_ticks.push_back(0);
                PPS_tick_correction.push_back(sumOfPreviousMissingTicks);
                continue;
            }
            
            if (missingPTicks != 0)
            {
                LAPPD_PPS_missing_ticks.push_back(missingPTicks);
                PPS_tick_correction.push_back(missingPTicks + sumOfPreviousMissingTicks);
                combined = true;
		std::cout << "Pushing PPS correction " << i << ", this PPS interval tick is " << thisInterval 
			<< ", missing ticks: " << thisMissingTicks << ", push missing " << LAPPD_PPS_missing_ticks.back() 
			<< ", push correction " << PPS_tick_correction.back() << std::endl;
                continue;
            }

            // Case 3: If one PPS is not recoreded correctly, like one interval is 1574262436, followed by a 4825737559
            // Then 1574262436 + 4825737559 = 6399999995 = 3200000000 * 2 - 5
            long long nextInterval = 0;

            if (i < LAPPD_PPS.size() - 1)
            {
                nextInterval = static_cast<long long>(LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE) - static_cast<long long>(thisPPS);
                long long combinedInterval = nextInterval + thisInterval;
		long long combinedRemainder = combinedInterval % static_cast<long long>(intervalTicks);

		// Case 3a: If the combined interval is an exact multiple of IntervalTicks
		if (combinedRemainder == 0) {
		    LAPPD_PPS_missing_ticks.push_back(0);
		    PPS_tick_correction.push_back(sumOfPreviousMissingTicks);
		    combined = true;

		    std::cout << "Pushing PPS correction " << i << ", this PPS interval tick is " << thisInterval 
			    << ", missing ticks: " << thisMissingTicks << ", push missing " << LAPPD_PPS_missing_ticks.back()
			    << ", push correction " << PPS_tick_correction.back() << std::endl;
		    continue;
		}

		// Case 3b: If the combined Interval is close to the IntervalTicks
		long long combinedMissingTicks = static_cast<long long>(intervalTicks) - combinedRemainder;
                if (combinedMissingTicks > -30 && combinedMissingTicks < 30)
                {
                    LAPPD_PPS_missing_ticks.push_back(combinedMissingTicks);
                    PPS_tick_correction.push_back(combinedMissingTicks + sumOfPreviousMissingTicks);
                    combined = true;
		    std::cout << "Pushing PPS correction " << i << ", this PPS interval tick is " << thisInterval 
			    << ", missing ticks: " << thisMissingTicks << ", push missing " << LAPPD_PPS_missing_ticks.back()
			    << ", push correction " << PPS_tick_correction.back() << std::endl;
                    continue;
                }
            }

            if (!combined)
            {
                LAPPD_PPS_missing_ticks.push_back(0);
                PPS_tick_correction.push_back(sumOfPreviousMissingTicks); 
	    }
        }   
	std::cout << "Pushing PPS correction " << i << ", this PPS interval tick is " << thisInterval 
		<< ", missing ticks: " << thisMissingTicks << ", push missing " << LAPPD_PPS_missing_ticks.back() 
		<< ", push correction " << PPS_tick_correction.back() << std::endl;
    }

    // loop all data events, plus the offset, save the event time and beamgate time
    // fing the closest CTC trigger, also save all information
    std::cout << "Start saving results. PPS size: " << LAPPD_PPS.size() << ", beamgate size: " << LAPPDDataBeamgateUL.size() << std::endl;
    std::cout << "First PPS: " << LAPPD_PPS.at(0) / LAPPD_PPS_SCALE << ", Last PPS: " << LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE << std::endl;
    std::cout << "First BG: " << LAPPDDataBeamgateUL.at(0) / LAPPD_PPS_SCALE 
	<< ", Last BG: " << LAPPDDataBeamgateUL.at(LAPPDDataBeamgateUL.size() - 1) / LAPPD_PPS_SCALE << std::endl;
    std::cout << "First TS: " << LAPPDDataTimeStampUL.at(0) / LAPPD_PPS_SCALE 
	<< ", Last TS: " << LAPPDDataTimeStampUL.at(LAPPDDataTimeStampUL.size() - 1) / LAPPD_PPS_SCALE << std::endl;

    for (int l = 0; l < LAPPDDataTimeStampUL.size(); l++)
    {
        ULong64_t TS_ns = LAPPDDataTimeStampUL.at(l) / 1000;
        ULong64_t TS_truncated_ps = LAPPDDataTimeStampUL.at(l) % 1000;
        ULong64_t BG_ns = LAPPDDataBeamgateUL.at(l) / 1000;
        ULong64_t BG_truncated_ps = LAPPDDataBeamgateUL.at(l) % 1000;
        ULong64_t driftCorrectionForTS = 0;
        ULong64_t driftCorrectionForBG = 0;

        if (drift != 0)
        {
            double driftScaling = static_cast<double>(TS_ns) / (trueInterval * 1.0 / 1000);
            driftCorrectionForTS = static_cast<ULong64_t>(static_cast<double>(drift) * driftScaling / 1000);
            double driftScalingBG = static_cast<double>(BG_ns) / (trueInterval * 1.0 / 1000);
            driftCorrectionForBG = static_cast<ULong64_t>(static_cast<double>(drift) * driftScalingBG / 1000);
            
	    std::cout << "drift = " << drift << ", driftScaling = " << driftScaling << ", driftCorrectionForTS = " << driftCorrectionForTS 
		<< ", driftCorrectionForBG = " << driftCorrectionForBG << std::endl;
            std::cout << "driftCorrectionForTS = " << driftCorrectionForTS << ", driftCorrectionForBG = " << driftCorrectionForBG << std::endl;
        }

        ULong64_t DriftCorrectedTS_ns = TS_ns + final_offset_ns + driftCorrectionForTS;
        ULong64_t DriftCorrectedBG_ns = BG_ns + final_offset_ns + driftCorrectionForBG;
        std::cout << "DriftCorrectedTS_ns = " << DriftCorrectedTS_ns << " =: final_offset_ns = " << final_offset_ns
	    << " + driftCorrectionForTS = " << driftCorrectionForTS << " + TS_ns = " << TS_ns << std::endl;

        ULong64_t min_mean_dev_match = std::numeric_limits<ULong64_t>::max();
        int matchedIndex = 0;
        for (int c = 0; c < CTCTrigger.size(); c++)
        {
            ULong64_t CTCTrigger_ns = CTCTrigger.at(c);
            Long64_t diff = static_cast<Long64_t>(CTCTrigger_ns) - static_cast<Long64_t>(DriftCorrectedTS_ns);
	    
	    if (diff < 0)
                diff = -diff;
            
	    if (static_cast<ULong64_t>(diff) < min_mean_dev_match)
            {
                matchedIndex = c;
                min_mean_dev_match = static_cast<ULong64_t>(diff);
            }
        }

        // check the position of beamgate and timestamp raw fall into which pps interval;
        bool TSFound = false;
        bool BGFound = false;
	const long long PPS_TICK_OFFSET = 1000;

        for (int i = 0; i < LAPPD_PPS.size() - 1; i++)
        {
            if (LAPPD_PPS.at(i) / LAPPD_PPS_SCALE < LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE && 
		LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE > LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE)
            {
                TimeStamp_correction_tick.push_back(PPS_tick_correction.at(i) + PPS_TICK_OFFSET);
                TSFound = true;
                break;
            }
        }
        
	if (!TSFound && LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE > LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE)
        {
            TimeStamp_correction_tick.push_back(PPS_tick_correction.at(LAPPD_PPS.size() - 1) + PPS_TICK_OFFSET);
            TSFound = true;
        }

        if (!TSFound && LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE < LAPPD_PPS.at(0) / LAPPD_PPS_SCALE)
        {
            TimeStamp_correction_tick.push_back(0 + PPS_TICK_OFFSET);
            TSFound = true;
        }

        for (int i = 0; i < LAPPD_PPS.size() - 1; i++)
        {
            if (LAPPD_PPS.at(i) / LAPPD_PPS_SCALE < LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE && 
		LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE > LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE)
            {
                BeamGate_correction_tick.push_back(PPS_tick_correction.at(i) + PPS_TICK_OFFSET);
                BGFound = true;
                std::cout << "Normal push: BG_raw = " << LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE << ", pps = " << LAPPD_PPS.at(i) 
			<< ", pps/" << LAPPD_PPS_SCALE << " = " << LAPPD_PPS.at(i) / LAPPD_PPS_SCALE << std::endl;
                if (LAPPD_PPS_interval_ticks.at(i) != intervalTicks)
                {
                    std::cout << "Warning: PPS interval is not " << intervalTicks << " at index " << i << ", it is " << LAPPD_PPS_interval_ticks.at(i) << std::endl;
                }
                break;
            }
        }
        
	if (!BGFound && LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE > LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE)
        {
            BeamGate_correction_tick.push_back(PPS_tick_correction.at(LAPPD_PPS.size() - 1) + PPS_TICK_OFFSET);
            std::cout << "BG_raw = " << LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE << ", pps = " << LAPPD_PPS.at(LAPPD_PPS.size() - 1) 
		<< ", pps/" << LAPPD_PPS_SCALE << " = " << LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE << std::endl;
            BGFound = true;
        }

        if (!BGFound && LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE < LAPPD_PPS.at(0) / LAPPD_PPS_SCALE)
        {
            BeamGate_correction_tick.push_back(0 + PPS_TICK_OFFSET);
            std::cout << "BG_raw less than pps0" << std::endl;
            BGFound = true;
        }

        if (!TSFound || !BGFound)
        {
            std::cout << "Error: PPS not found for event " << l << ", TSFound: " << TSFound << ", BGFound: " << BGFound << std::endl;
        }

        /*
        cout << "******Found result:" << endl;
        cout << "LAPPDDataTimeStampUL.at(" << l << "): " << LAPPDDataTimeStampUL.at(l)/3125 << endl;
        cout << "LAPPDDataBeamgateUL.at(" << l << "): " << LAPPDDataBeamgateUL.at(l)/3125 << endl;
        cout << "TS_ns: " << TS_ns << endl;
        cout << "BG_ns: " << BG_ns << endl;
        cout << "final_offset_ns: " << final_offset_ns << endl;
        cout << "drift: " << drift << endl;
        cout << "TS driftscaling: " << TS_ns / trueInterval / 1000 << endl;
        cout << "DriftCorrectedTS_ns: " << DriftCorrectedTS_ns << endl;
        cout << "TS_truncated_ps: " << TS_truncated_ps << endl;
        cout << "DriftCorrectedBG_ns: " << DriftCorrectedBG_ns << endl;
        cout << "BG_truncated_ps: " << BG_truncated_ps << endl;
        cout << "Found min_mean_dev_match: " << min_mean_dev_match << endl;
        cout << "MatchedIndex: " << matchedIndex << endl;
        cout << "CTCTrigger.at(" << matchedIndex << "): " << CTCTrigger.at(matchedIndex) << endl;
        cout << "BeamGate_correction_tick at " << l << ": " << BeamGate_correction_tick.at(l) << endl;
        cout << "TimeStamp_correction_tick at " << l << ": " << TimeStamp_correction_tick.at(l) << endl;
        */

        // for this LAPPDDataBeamgateUL.at(l), in the LAPPD_PPS vector, find it's closest PPS before and after, 
	// and also calculate the time difference between them, and the missing tick between them.
        // save as BG_PPSBefore_tick, BG_PPSAfter_tick, BG_PPSDiff_tick, BG_PPSMissing_tick
        // Do the same thing for LAPPDDataTimeStampUL.at(l), save as TS_PPSBefore_tick, TS_PPSAfter_tick, TS_PPSDiff_tick, TS_PPSMissing_tick

        ULong64_t BG_PPSBefore_tick = 0;
        ULong64_t BG_PPSAfter_tick = 0;
        ULong64_t BG_PPSDiff_tick = 0;
        ULong64_t BG_PPSMissing_tick = 0;

        ULong64_t TS_PPSBefore_tick = 0;
        ULong64_t TS_PPSAfter_tick = 0;
        ULong64_t TS_PPSDiff_tick = 0;
        ULong64_t TS_PPSMissing_tick = 0;
        
	// if the first PPS is later than beamgate, then set before = 0, after is the first, diff is the first, missing is the first
        std::cout << "Start finding PPS before and after the beamgate" << std::endl;

        if ((LAPPD_PPS.at(0) / LAPPD_PPS_SCALE > LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE) || 
		(LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE < LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE))
        {
            if (LAPPD_PPS.at(0) / LAPPD_PPS_SCALE > LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE)
            {
                BG_PPSBefore_tick = 0;
                BG_PPSAfter_tick = LAPPD_PPS.at(0) / LAPPD_PPS_SCALE;
                BG_PPSDiff_tick = LAPPD_PPS.at(0) / LAPPD_PPS_SCALE;
                BG_PPSMissing_tick = LAPPD_PPS.at(0) / LAPPD_PPS_SCALE;

                std::cout << "First PPS is later than beamgate, before is 0, after is the first, diff is the first, missing is the first" << std::endl;
            }

            // if the last PPS is earlier than beamgate, then set before is the last, after = 0, diff is the last, missing is the last
            if (LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE < LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE)
            {
                BG_PPSBefore_tick = LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE;
                BG_PPSAfter_tick = 0;
                BG_PPSDiff_tick = LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE;
                BG_PPSMissing_tick = LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE;
                
		std::cout << "Last PPS is earlier than beamgate, before is the last, after is 0, diff is the last, missing is the last" << std::endl;
            }
        }

        else
        {
            // if the first PPS is earlier than beamgate, and the last PPS is later than beamgate, then find the closest PPS before and after
            for (int i = 0; i < LAPPD_PPS.size() - 1; i++)
            {
                if (LAPPD_PPS.at(i) / LAPPD_PPS_SCALE < LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE && 
			LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE > LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE)
                {
                    BG_PPSBefore_tick = LAPPD_PPS.at(i) / LAPPD_PPS_SCALE;
                    BG_PPSAfter_tick = LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE;
                    BG_PPSDiff_tick = LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE - LAPPD_PPS.at(i) / LAPPD_PPS_SCALE;
                    ULong64_t DiffTick = (LAPPD_PPS.at(i + 1) - LAPPD_PPS.at(i)) / LAPPD_PPS_SCALE - PPS_TICK_OFFSET;

		    ULong64_t remainder = DiffTick % intervalTicks;
		    
		    // Perfect multiple of intervalTicks
		    if (remainder == 0) {
			BG_PPSMissing_tick = 0;
		    }
		    else {
			BG_PPSMissing_tick = std::min(remainder, intervalTicks - remainder);
 		    }

                    std::cout << "Found PPS before and after the beamgate, before: " << BG_PPSBefore_tick << ", after: " << BG_PPSAfter_tick 
			<< ", diff: " << BG_PPSDiff_tick << ", missing: " << BG_PPSMissing_tick << std::endl;
                    break;
                }
            }
        }

        // do the same for timestamp
        std::cout << "Start finding PPS before and after the timestamp" << std::endl;

        if ((LAPPD_PPS.at(0) / LAPPD_PPS_SCALE > LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE) || 
		(LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE < LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE))
        {
            if (LAPPD_PPS.at(0) / LAPPD_PPS_SCALE > LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE)
            {
                TS_PPSBefore_tick = 0;
                TS_PPSAfter_tick = LAPPD_PPS.at(0) / LAPPD_PPS_SCALE;
                TS_PPSDiff_tick = LAPPD_PPS.at(0) / LAPPD_PPS_SCALE;
                TS_PPSMissing_tick = LAPPD_PPS.at(0) / LAPPD_PPS_SCALE;

                std::cout << "First PPS is later than timestamp, before is 0, after is the first, diff is the first, missing is the first" << std::endl;
            }

            if (LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE < LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE)
            {
                TS_PPSBefore_tick = LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE;
                TS_PPSAfter_tick = 0;
                TS_PPSDiff_tick = LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE;
                TS_PPSMissing_tick = LAPPD_PPS.at(LAPPD_PPS.size() - 1) / LAPPD_PPS_SCALE;

                std::cout << "Last PPS is earlier than timestamp, before is the last, after is 0, diff is the last, missing is the last" << std::endl;
            }
        }

        else
        {
            for (int i = 0; i < LAPPD_PPS.size() - 1; i++)
            {
                if (LAPPD_PPS.at(i) / LAPPD_PPS_SCALE < LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE && 
			LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE > LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE)
                {
                    TS_PPSBefore_tick = LAPPD_PPS.at(i) / LAPPD_PPS_SCALE;
                    TS_PPSAfter_tick = LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE;
                    TS_PPSDiff_tick = LAPPD_PPS.at(i + 1) / LAPPD_PPS_SCALE - LAPPD_PPS.at(i) / LAPPD_PPS_SCALE;
                    ULong64_t DiffTick = (LAPPD_PPS.at(i + 1) - LAPPD_PPS.at(i)) / LAPPD_PPS_SCALE - PPS_TICK_OFFSET;

		    ULong64_t remainder = DiffTick % intervalTicks;

                    // Perfect multiple of intervalTicks
                    if (remainder == 0) {
                        TS_PPSMissing_tick = 0;
                    }
                    else {
                        TS_PPSMissing_tick = std::min(remainder, intervalTicks - remainder);
                    }

                    std::cout << "Found PPS before and after the timestamp, before: " << TS_PPSBefore_tick << ", after: " << TS_PPSAfter_tick 
			<< ", diff: " << TS_PPSDiff_tick << ", missing: " << TS_PPSMissing_tick << std::endl;
                    break;
                }
            }
        }

        TimeStampRaw.push_back(LAPPDDataTimeStampUL.at(l) / LAPPD_PPS_SCALE);
        BeamGateRaw.push_back(LAPPDDataBeamgateUL.at(l) / LAPPD_PPS_SCALE);
        TimeStamp_ns.push_back(DriftCorrectedTS_ns);
        BeamGate_ns.push_back(DriftCorrectedBG_ns);
        TimeStamp_ps.push_back(TS_truncated_ps);
        BeamGate_ps.push_back(BG_truncated_ps);
        EventIndex.push_back(l);
        EventDeviation_ns.push_back(min_mean_dev_match);
        // to get ps, should minus the TS_truncated_ps
        CTCTriggerIndex.push_back(matchedIndex);
        CTCTriggerTimeStamp_ns.push_back(CTCTrigger.at(matchedIndex));

        BG_PPSBefore.push_back(BG_PPSBefore_tick);
        BG_PPSAfter.push_back(BG_PPSAfter_tick);
        BG_PPSDiff.push_back(BG_PPSDiff_tick);
        BG_PPSMiss.push_back(BG_PPSMissing_tick);
        TS_PPSBefore.push_back(TS_PPSBefore_tick);
        TS_PPSAfter.push_back(TS_PPSAfter_tick);
        TS_PPSDiff.push_back(TS_PPSDiff_tick);
        TS_PPSMiss.push_back(TS_PPSMissing_tick);

        TS_driftCorrection_ns.push_back(driftCorrectionForTS);
        BG_driftCorrection_ns.push_back(driftCorrectionForBG);
    }

    ULong64_t totalEventNumber = LAPPDDataTimeStampUL.size();
    ULong64_t gotOrphanCount_out = gotOrphanCount;
    ULong64_t gotMin_mean_dev_noOrphan_out = gotMin_mean_dev_noOrphan;
    ULong64_t increament_times_out = increament_times;
    ULong64_t min_mean_dev_out = min_mean_dev;
    ULong64_t final_i_out = final_i;
    ULong64_t final_j_out = final_j;
    ULong64_t drift_out = drift;

    vector<ULong64_t> FitInfo = {final_offset_ns, final_offset_ps_negative, gotOrphanCount_out, gotMin_mean_dev_noOrphan_out, 
	increament_times_out, min_mean_dev_out, final_i_out, final_j_out, totalEventNumber, drift_out};

    vector<vector<ULong64_t>> Result = {FitInfo, TimeStampRaw, BeamGateRaw, TimeStamp_ns, BeamGate_ns, TimeStamp_ps, BeamGate_ps, 
	EventIndex, EventDeviation_ns, CTCTriggerIndex, CTCTriggerTimeStamp_ns, BeamGate_correction_tick, TimeStamp_correction_tick, 
	LAPPD_PPS_interval_ticks, BG_PPSBefore, BG_PPSAfter, BG_PPSDiff, BG_PPSMiss, TS_PPSBefore, TS_PPSAfter, TS_PPSDiff, TS_PPSMiss, 
	TS_driftCorrection_ns, BG_driftCorrection_ns};
   
    return Result;
}

vector<vector<ULong64_t>> fitInPartFile(TTree *lappdTree, TTree *triggerTree, int runNumber, int subRunNumber, int partFileNumber, 
		int LAPPD_ID, int fitTargetTriggerWord, bool triggerGrouped, int intervalInSecond)
{
    std::cout << "***************************************" << std::endl;
    std::cout << "Fitting in run " << runNumber << ", sub run " << subRunNumber << ", part file " << partFileNumber << " for LAPPD ID " << LAPPD_ID << std::endl;

    vector<ULong64_t> LAPPD_PPS0;
    vector<ULong64_t> LAPPD_PPS1;
    vector<ULong64_t> LAPPDDataTimeStamp0_UL;
    vector<ULong64_t> LAPPDDataBeamgate0_UL;
    vector<ULong64_t> LAPPDDataTimeStamp1_UL;
    vector<ULong64_t> LAPPDDataBeamgate1_UL;

    vector<ULong64_t> CTCTargetTimeStamp;
    vector<ULong64_t> CTCPPSTimeStamp;

    int LAPPD_ID_inTree;
    int runNumber_inTree;
    int subRunNumber_inTree;
    int partFileNumber_inTree;

    ULong64_t ppsTime0;
    ULong64_t ppsTime1;
    ULong64_t LAPPDTimeStamp0_UL;
    ULong64_t LAPPDBeamgate0_UL;
    ULong64_t LAPPDTimeStamp1_UL;
    ULong64_t LAPPDBeamgate1_UL;

    vector<uint32_t> *CTCTriggerWord = nullptr;
    ULong64_t CTCTimeStamp;
    vector<uint32_t> *groupedTriggerWords = nullptr;
    vector<ULong64_t> *groupedTriggerTimestamps = nullptr;

    lappdTree->SetBranchAddress("RunNumber", &runNumber_inTree);
    lappdTree->SetBranchAddress("SubRunNumber", &subRunNumber_inTree);
    lappdTree->SetBranchAddress("PartFileNumber", &partFileNumber_inTree);
    lappdTree->SetBranchAddress("LAPPD_ID", &LAPPD_ID_inTree);
    lappdTree->SetBranchAddress("LAPPDDataTimeStamp0_UL", &LAPPDTimeStamp0_UL);
    lappdTree->SetBranchAddress("LAPPDDataBeamgate0_UL", &LAPPDBeamgate0_UL);
    lappdTree->SetBranchAddress("LAPPDDataTimeStamp1_UL", &LAPPDTimeStamp1_UL);
    lappdTree->SetBranchAddress("LAPPDDataBeamgate1_UL", &LAPPDBeamgate1_UL);
    lappdTree->SetBranchAddress("ppsTime0", &ppsTime0);
    lappdTree->SetBranchAddress("ppsTime1", &ppsTime1);

    // triggerTree->Print();
    triggerTree->SetBranchAddress("RunNumber", &runNumber_inTree);
    triggerTree->SetBranchAddress("SubRunNumber", &subRunNumber_inTree);
    triggerTree->SetBranchAddress("PartFileNumber", &partFileNumber_inTree);
    
    if (triggerGrouped)
    {
        triggerTree->SetBranchAddress("gTrigWord", &groupedTriggerWords);
        triggerTree->SetBranchAddress("gTrigTime", &groupedTriggerTimestamps);
    }
    else
    {
        triggerTree->SetBranchAddress("CTCTriggerWord", &CTCTriggerWord);
        triggerTree->SetBranchAddress("CTCTimeStamp", &CTCTimeStamp);
    }

    int l_nEntries = lappdTree->GetEntries();
    int t_nEntries = triggerTree->GetEntries();
    int repeatedPPSNumber0 = -1;
    int repeatedPPSNumber1 = -1;

    // 4. For each LAPPD ID,
    // 4a. Loop over TimeStamp tree
    for (int i = 0; i < l_nEntries; i++)
    {
        lappdTree->GetEntry(i);
        
	// 1000 / 8 * 25 = 3125
        if (LAPPD_ID_inTree == LAPPD_ID && runNumber_inTree == runNumber && subRunNumber_inTree == subRunNumber && partFileNumber_inTree == partFileNumber)
        {
            // Define boolean flags for readability
            bool hasData0 = (LAPPDTimeStamp0_UL != 0);
            bool hasData1 = (LAPPDTimeStamp1_UL != 0);

            // --- 1. Handle Data Events (Strict Coincidence) ---
            if (hasData0 && hasData1)
            {
                // Both boards have valid data. Safe to push!
                
                std::cout << "In unit of ps, after conversion and saving, LAPPDTimeStamp0_UL: " << LAPPDTimeStamp0_UL * 3125 
			<< ", LAPPDBeamgate0_UL: " << LAPPDBeamgate0_UL * 3125 << std::endl;
                std::cout << "In second, use double, Timestamp0: " << static_cast<double>(LAPPDTimeStamp0_UL * 3125) / 1E12 
			<< ", Beamgate0: " << static_cast<double>(LAPPDBeamgate0_UL * 3125) / 1E12 << std::endl;
               
                 
                std::cout << "In unit of ps, after conversion and saving, LAPPDTimeStamp1_UL: " << LAPPDTimeStamp1_UL * 3125 
			<< ", LAPPDBeamgate1_UL: " << LAPPDBeamgate1_UL * 3125 << std::endl;
                std::cout << "In second, use double, Timestamp1: " << static_cast<double>(LAPPDTimeStamp1_UL * 3125) / 1E12 
			<< ", Beamgate1: " << static_cast<double>(LAPPDBeamgate1_UL * 3125) / 1E12 << std::endl;
                
                LAPPDDataTimeStamp0_UL.push_back(LAPPDTimeStamp0_UL * 3125);
                LAPPDDataBeamgate0_UL.push_back(LAPPDBeamgate0_UL * 3125);
                
                LAPPDDataTimeStamp1_UL.push_back(LAPPDTimeStamp1_UL * 3125);
                LAPPDDataBeamgate1_UL.push_back(LAPPDBeamgate1_UL * 3125);
            }
            else if (hasData0 || hasData1)
            {
                // Mismatch: One board fired, the other did not.
                // We drop the event to keep the vectors perfectly parallel.
                std::cout << "Warning: Mismatched data at tree index " << i 
                          << ". Dropping event to maintain ACDC alignment." << std::endl;
            }

            // --- 2. Handle PPS Heartbeats Independently ---
            if (!hasData0)
            {
                if (LAPPD_PPS0.size() == 0)
                    LAPPD_PPS0.push_back(ppsTime0 * 3125);
                else if (ppsTime0 * 3125 != LAPPD_PPS0.at(LAPPD_PPS0.size() - 1))
                    LAPPD_PPS0.push_back(ppsTime0 * 3125);
                else
                    repeatedPPSNumber0 += 1;
            }
            
            if (!hasData1)
            {
                if (LAPPD_PPS1.size() == 0)
                    LAPPD_PPS1.push_back(ppsTime1 * 3125);    
                else if (ppsTime1 * 3125 != LAPPD_PPS1.at(LAPPD_PPS1.size() - 1))
                    LAPPD_PPS1.push_back(ppsTime1 * 3125);
                else
                    repeatedPPSNumber1 += 1;
            }
        }
    }
    
    std::cout << "repeatedPPSNumber0: " << repeatedPPSNumber0 << ", loaded PPS0 size: " << LAPPD_PPS0.size() << std::endl;
    std::cout << "repeatedPPSNumber1: " << repeatedPPSNumber1 << ", loaded PPS1 size: " << LAPPD_PPS1.size() << std::endl;

    // 4b. Loop over Trig (or GTrig) tree
    for (int i = 0; i < t_nEntries; i++)
    {
        triggerTree->GetEntry(i);
        if (runNumber_inTree == runNumber && subRunNumber_inTree == subRunNumber && partFileNumber_inTree == partFileNumber)
        {
            // std::cout << "triggerGrouped: " << triggerGrouped << ", fitTargetTriggerWord: " << fitTargetTriggerWord << std::endl;
            if (triggerGrouped)
            {
                for (int j = 0; j < groupedTriggerWords->size(); j++)
                {
                    // std::cout << "At j = " << j << ",finding groupedTriggerWords: " << groupedTriggerWords->at(j) 
			// << ", fitTargetTriggerWord: " << fitTargetTriggerWord << std::endl;

                    if (groupedTriggerWords->at(j) == fitTargetTriggerWord)
                        CTCTargetTimeStamp.push_back(groupedTriggerTimestamps->at(j));
                    if (groupedTriggerWords->at(j) == 32)
                        CTCPPSTimeStamp.push_back(groupedTriggerTimestamps->at(j));
                }
            }
            else
            {
                for (int j = 0; j < CTCTriggerWord->size(); j++)
                {
                    if (CTCTriggerWord->at(j) == fitTargetTriggerWord)
                        CTCTargetTimeStamp.push_back(CTCTimeStamp);
                    if (CTCTriggerWord->at(j) == 32)
                        CTCPPSTimeStamp.push_back(CTCTimeStamp);
                }
            }
        }
    }
    
    std::cout << "Vector for partfile " << partFileNumber << " for LAPPD ID " << LAPPD_ID << " loaded." << std::endl;
    std::cout << "LAPPDDataTimeStamp0_UL in ps size: " << LAPPDDataTimeStamp0_UL.size() << std::endl;
    std::cout << "LAPPDDataTimeStamp1_UL in ps size: " << LAPPDDataTimeStamp1_UL.size() << std::endl;
    std::cout << "LAPPDDataBeamgate0_UL in ps size: " << LAPPDDataBeamgate0_UL.size() << std::endl;
    std::cout << "LAPPDDataBeamgate1_UL in ps size: " << LAPPDDataBeamgate1_UL.size() << std::endl;
    std::cout << "LAPPD_PPS0 size: " << LAPPD_PPS0.size() << std::endl;
    std::cout << "LAPPD_PPS1 size: " << LAPPD_PPS1.size() << std::endl;
    std::cout << "CTCTargetTimeStamp size: " << CTCTargetTimeStamp.size() << std::endl;
    std::cout << "CTCPPSTimeStamp size: " << CTCPPSTimeStamp.size() << std::endl;
   
    // 5. Find the number of resets in LAPPD PPS: 
    
    // Initialize the fit stop index to the last valid timestamp index.
    // If no reset is found, the full timestamp vector will be used.
    int LAPPDDataFitStopIndex0 = LAPPDDataTimeStamp0_UL.size() - 1;
    int LAPPDDataFitStopIndex1 = LAPPDDataTimeStamp1_UL.size() - 1;
    int resetNumber0 = 0;
    int resetNumber1 = 0;
  
    // Check ACDC 0 and ACDC 1 independently for a PPS reset. 
    // If all PPS in this part file are increasing, then there is no reset
    for (int i = 1; i < LAPPD_PPS0.size(); i++)
    {
        if (LAPPD_PPS0[i] < LAPPD_PPS0[i - 1])
        {
            resetNumber0 += 1;
            std::cout << "For LAPPD ID " << LAPPD_ID << ", run number " << runNumber << ", sub run number " << subRunNumber 
		<< ", part file number " << partFileNumber << ", reset " << " found at PPS_ACDC0 index " << i << std::endl;
            break;
        }
    }
    for (int i = 1; i < LAPPD_PPS1.size(); i++)
    {
        if (LAPPD_PPS1[i] < LAPPD_PPS1[i - 1])
        {
            resetNumber1 += 1;
            std::cout << "For LAPPD ID " << LAPPD_ID << ", run number " << runNumber << ", sub run number " << subRunNumber 
		<< ", part file number " << partFileNumber << ", reset " << " found at PPS_ACDC1 index " << i << std::endl;
            break;
        }
    }

    if (resetNumber0 == 0) {
        std::cout << "For LAPPD ID " << LAPPD_ID << ", run number " << runNumber << ", sub run number " << subRunNumber 
		<< ", part file number " << partFileNumber << ", no reset found for ACDC 0." << std::endl;
    }
    if (resetNumber1 == 0) {
        std::cout << "For LAPPD ID " << LAPPD_ID << ", run number " << runNumber << ", sub run number " << subRunNumber 
		<< ", part file number " << partFileNumber << ", no reset found for ACDC 1." << std::endl;
    }

    // If a reset is found, determine the fit stop index separately for each ACDC.
    if (resetNumber0 != 0)
    {
        for (int i = 1; i < LAPPDDataTimeStamp0_UL.size(); i++)
        {
            if (LAPPDDataTimeStamp0_UL[i] < LAPPDDataTimeStamp0_UL[i - 1])
            {
                LAPPDDataFitStopIndex0 = i - 1;
                break;
            }
        }
        // TODO: extend this to later offsets
    }
    if (resetNumber1 != 0)
    {
        for (int i = 1; i < LAPPDDataTimeStamp1_UL.size(); i++)
        {
            if (LAPPDDataTimeStamp1_UL[i] < LAPPDDataTimeStamp1_UL[i - 1])
            {
                LAPPDDataFitStopIndex1 = i - 1;
                break;
            }
        }
        // TODO: extend this to later offsets
    }
    
    // 6. Use the target trigger word to fit the offset
    // TODO: Fit for each reset
    // Note by Anuj Gupta: I have checked all of the PPS Timestamps for runs <= 5668 and no reset has been observed.

    vector<vector<ULong64_t>> ResultTotal;

    // Sanity checks for timing vector sizes across ACDCs
    if (LAPPDDataTimeStamp0_UL.size() != LAPPDDataTimeStamp1_UL.size()) {
        std::cout << "Warning: ACDC 0 and ACDC 1 timestamp vector sizes differ: " << LAPPDDataTimeStamp0_UL.size() << " vs "<< LAPPDDataTimeStamp1_UL.size() << std::endl;
    }
    if (LAPPDDataBeamgate0_UL.size() != LAPPDDataBeamgate1_UL.size()) {
        std::cout << "Warning: ACDC 0 and ACDC 1 beamgate vector sizes differ: " << LAPPDDataBeamgate0_UL.size() << " vs "<< LAPPDDataBeamgate1_UL.size() << std::endl;
    }
    if (LAPPD_PPS0.size() != LAPPD_PPS1.size()) {
        std::cout << "Warning: ACDC 0 and ACDC 1 PPS vector sizes differ: " << LAPPD_PPS0.size() << " vs " << LAPPD_PPS1.size() << std::endl;
    }

    if (LAPPDDataTimeStamp0_UL.size() == LAPPDDataBeamgate0_UL.size() && 
        LAPPDDataTimeStamp1_UL.size() == LAPPDDataBeamgate1_UL.size())
    {
        if (LAPPD_PPS0.size() == 0 || LAPPD_PPS1.size() == 0)
        {
            std::cout << "Error: PPS0 or PPS1 is empty, return empty result." << std::endl;
            return ResultTotal;
        }

        vector<vector<ULong64_t>> ResultACDC0 = fitInThisReset(LAPPDDataTimeStamp0_UL, LAPPDDataBeamgate0_UL, LAPPD_PPS0, fitTargetTriggerWord, 
		CTCTargetTimeStamp, CTCPPSTimeStamp, static_cast<ULong64_t>(intervalInSecond) * 1000000000000ULL, partFileNumber, LAPPD_ID, 0);
        vector<vector<ULong64_t>> ResultACDC1 = fitInThisReset(LAPPDDataTimeStamp1_UL, LAPPDDataBeamgate1_UL, LAPPD_PPS1, fitTargetTriggerWord, 
		CTCTargetTimeStamp, CTCPPSTimeStamp, static_cast<ULong64_t>(intervalInSecond) * 1000000000000ULL, partFileNumber, LAPPD_ID, 1);

        // 7. Save the offset for this LAPPD ID, run number, part file number, index and reset number.
        std::cout << "Fitting in part file " << partFileNumber << " for LAPPD ID " << LAPPD_ID << " done." << std::endl;

        // Combine ResultACDC0 and ResultACDC1 to ResultTotal
        for (int i = 0; i < ResultACDC0.size(); i++)
        {
            ResultTotal.push_back(ResultACDC0[i]);
        }
        for (int i = 0; i < ResultACDC1.size(); i++)
        {
            ResultTotal.push_back(ResultACDC1[i]);
        }
    }
    return ResultTotal;
}

void offsetFit_MultipleLAPPD(string fileName, int fitTargetTriggerWord, bool triggerGrouped, int intervalInSecond, int processPFNumber)
{

    TH1::AddDirectory(kFALSE);

    // 1. Load LAPPDTree.root
    const string file = fileName;
    TFile *f = new TFile(file.c_str(), "READ");
    
    if (!f->IsOpen())
    {
        std::cerr << "Error: cannot open file " << file << std::endl;
        return;
    }
    std::cout << "Opened file " << file << std::endl;

    std::ofstream outputOffset("offset.txt");
    outputOffset << "runNumber"
                 << "\t"
                 << "subRunNumber"
                 << "\t"
                 << "partFileNumber"
                 << "\t"
                 << "resetNumber"
                 << "\t"
                 << "LAPPD_ID"
                 << "\t"
                 << "offset_ACDC0_ns"
                 << "\t"
                 << "offset_ACDC1_ns"
                 << "\t"
                 << "offset_ACDC0_ps_negative"
                 << "\t"
                 << "offset_ACDC1_ps_negative"
                 << "\t"
                 << "gotOrphanCount_ACDC0"
                 << "\t"
                 << "gotOrphanCount_ACDC1"
                 << "\t"
                 << "n_beamgates"
                 << "\t"
                 << "min_mean_dev_noOrphan_ACDC0"
                 << "\t"
                 << "min_mean_dev_noOrphan_ACDC1"
                 << "\t"
                 << "increment_times_ACDC0"
                 << "\t"
                 << "increment_times_ACDC1"
                 << "\t"
                 << "min_mean_dev_ACDC0"
                 << "\t"
                 << "min_mean_dev_ACDC1"
                 << std::endl;

    // 2. Load TimeStamp tree, get run number and part file number
    TTree *TSTree = (TTree *)f->Get("TimeStamp");
    TTree *CTCTree = (TTree *)f->Get("Trig");
    TTree *GCTCTree = (TTree *)f->Get("GTrig");

    int runNumber;
    int subRunNumber;
    int partFileNumber;
    int LAPPD_ID;
    ULong64_t ppsCount0;
    
    TSTree->SetBranchAddress("RunNumber", &runNumber);
    TSTree->SetBranchAddress("SubRunNumber", &subRunNumber);
    TSTree->SetBranchAddress("PartFileNumber", &partFileNumber);
    TSTree->SetBranchAddress("LAPPD_ID", &LAPPD_ID);
    TSTree->SetBranchAddress("ppsCount0", &ppsCount0);
    
    std::vector<vector<int>> loopInfo;
    int nEntries = TSTree->GetEntries();
    
    for (int i = 0; i < nEntries; i++)
    {
        TSTree->GetEntry(i);
        
        // If this entry is PPS event, continue.
        if (ppsCount0 != 0)
            continue;
        
        vector<int> info = {runNumber, subRunNumber, partFileNumber, LAPPD_ID};
        
        // Find if this info is already in loopInfo
        bool found = std::find_if(loopInfo.begin(), loopInfo.end(), [&info](const std::vector<int> &vec)
                                  {
                                      return vec == info; // compare vec and info
                                  }) != loopInfo.end();   // if find_if doesn't return end，found

        if (!found) loopInfo.push_back(info);
    }

    // 3. For each unique run number and part file number, perform the fit.
    std::map<string, vector<vector<ULong64_t>>> ResultMap;
    std::vector<vector<int>>::iterator it;
    int pfNumber = 0;

    for (it = loopInfo.begin(); it != loopInfo.end(); it++)
    {
        int runNumber = (*it)[0];
        int subRunNumber = (*it)[1];
        int partFileNumber = (*it)[2];
        int LAPPD_ID = (*it)[3];
        
        if (processPFNumber != 0)
        {
            if (pfNumber >= processPFNumber)
                break;
        }

        vector<vector<ULong64_t>> Result;
        if (!triggerGrouped)
        { 
            // std::cout << "Trigger is not grouped" << std::endl;
            Result = fitInPartFile(TSTree, CTCTree, runNumber, subRunNumber, partFileNumber, LAPPD_ID, fitTargetTriggerWord, triggerGrouped, intervalInSecond);
        }
        else
        { 
            // std::cout << "Trigger is grouped" << std::endl;
            Result = fitInPartFile(TSTree, GCTCTree, runNumber, subRunNumber, partFileNumber, LAPPD_ID, fitTargetTriggerWord, triggerGrouped, intervalInSecond);
        }

        // Combine the *it to a string, and save the result to ResultMap
        string key = std::to_string(runNumber) + "_" + std::to_string(subRunNumber) + "_" + std::to_string(partFileNumber) + "_" + std::to_string(LAPPD_ID);
        ResultMap[key] = Result;
        pfNumber++;
    }

    // Loop the ResultMap, save the result to a root tree in a root file
    std::cout << "Start saving the result to root file and txt file..." << std::endl;

    TFile *fOut = new TFile("offsetFitResult.root", "RECREATE");
    TTree *tOut = new TTree("Events", "Events");
    
    int runNumber_out;
    int subRunNumber_out;
    int partFileNumber_out;
    int LAPPD_ID_out;
    ULong64_t EventIndex;
    ULong64_t EventNumberInThisPartFile; 
    ULong64_t final_offset_ns_0;
    ULong64_t final_offset_ns_1;
    ULong64_t final_offset_ps_negative_0;
    ULong64_t final_offset_ps_negative_1;
    double drift0;
    double drift1;
    ULong64_t gotOrphanCount_0;
    ULong64_t gotOrphanCount_1;
    ULong64_t gotMin_mean_dev_noOrphan_0;
    ULong64_t gotMin_mean_dev_noOrphan_1;
    ULong64_t increament_times_0;
    ULong64_t increament_times_1;
    ULong64_t min_mean_dev_0;
    ULong64_t min_mean_dev_1;
    ULong64_t TimeStampRaw_0, TimeStampRaw_1;
    ULong64_t BeamGateRaw_0, BeamGateRaw_1;
    ULong64_t TimeStamp_ns_0, TimeStamp_ns_1;
    ULong64_t BeamGate_ns_0, BeamGate_ns_1;
    ULong64_t TimeStamp_ps_0, TimeStamp_ps_1;
    ULong64_t BeamGate_ps_0, BeamGate_ps_1;
    ULong64_t EventDeviation_ns_0;
    ULong64_t EventDeviation_ns_1;
    ULong64_t CTCTriggerIndex_0, CTCTriggerIndex_1;
    ULong64_t CTCTriggerTimeStamp_ns_0, CTCTriggerTimeStamp_ns_1;
    long long BGMinusTrigger_ns_0, BGMinusTrigger_ns_1;
    long long BGCorrection_tick_0, BGCorrection_tick_1;
    long long TSCorrection_tick_0, TSCorrection_tick_1;
    ULong64_t LAPPD_PPS_interval_ticks_0, LAPPD_PPS_interval_ticks_1;
    ULong64_t BG_PPSBefore_tick_0, BG_PPSBefore_tick_1;
    ULong64_t BG_PPSAfter_tick_0, BG_PPSAfter_tick_1;
    ULong64_t BG_PPSDiff_tick_0, BG_PPSDiff_tick_1;
    ULong64_t BG_PPSMissing_tick_0, BG_PPSMissing_tick_1;
    ULong64_t TS_PPSBefore_tick_0, TS_PPSBefore_tick_1;
    ULong64_t TS_PPSAfter_tick_0, TS_PPSAfter_tick_1;
    ULong64_t TS_PPSDiff_tick_0, TS_PPSDiff_tick_1;
    ULong64_t TS_PPSMissing_tick_0, TS_PPSMissing_tick_1;
    ULong64_t TS_driftCorrection_ns_0, TS_driftCorrection_ns_1;
    ULong64_t BG_driftCorrection_ns_0, BG_driftCorrection_ns_1;

    tOut->Branch("runNumber", &runNumber_out, "runNumber/I");
    tOut->Branch("subRunNumber", &subRunNumber_out, "subRunNumber/I");
    tOut->Branch("partFileNumber", &partFileNumber_out, "partFileNumber/I");
    tOut->Branch("LAPPD_ID", &LAPPD_ID_out, "LAPPD_ID/I");
    tOut->Branch("EventIndex", &EventIndex, "EventIndex/l");
    tOut->Branch("n_beamgates", &EventNumberInThisPartFile, "n_beamgates/l");
    tOut->Branch("final_offset_ns_0", &final_offset_ns_0, "final_offset_ns_0/l");
    tOut->Branch("final_offset_ns_1", &final_offset_ns_1, "final_offset_ns_1/l");
    tOut->Branch("final_offset_ps_negative_0", &final_offset_ps_negative_0, "final_offset_ps_negative_0/l");
    tOut->Branch("final_offset_ps_negative_1", &final_offset_ps_negative_1, "final_offset_ps_negative_1/l");
    tOut->Branch("drift0", &drift0, "drift0/D");
    tOut->Branch("drift1", &drift1, "drift1/D");
    tOut->Branch("gotOrphanCount_0", &gotOrphanCount_0, "gotOrphanCount_0/l");
    tOut->Branch("gotOrphanCount_1", &gotOrphanCount_1, "gotOrphanCount_1/l");
    tOut->Branch("gotMin_mean_dev_noOrphan_0", &gotMin_mean_dev_noOrphan_0, "gotMin_mean_dev_noOrphan_0/l");
    tOut->Branch("gotMin_mean_dev_noOrphan_1", &gotMin_mean_dev_noOrphan_1, "gotMin_mean_dev_noOrphan_1/l");
    tOut->Branch("increament_times_0", &increament_times_0, "increament_times_0/l");
    tOut->Branch("increament_times_1", &increament_times_1, "increament_times_1/l");
    tOut->Branch("min_mean_dev_0", &min_mean_dev_0, "min_mean_dev_0/l");
    tOut->Branch("min_mean_dev_1", &min_mean_dev_1, "min_mean_dev_1/l");
    tOut->Branch("TimeStampRaw_0", &TimeStampRaw_0, "TimeStampRaw_0/l");
    tOut->Branch("TimeStampRaw_1", &TimeStampRaw_1, "TimeStampRaw_1/l");
    tOut->Branch("BeamGateRaw_0", &BeamGateRaw_0, "BeamGateRaw_0/l");
    tOut->Branch("BeamGateRaw_1", &BeamGateRaw_1, "BeamGateRaw_1/l");
    tOut->Branch("TimeStamp_ns_0", &TimeStamp_ns_0, "TimeStamp_ns_0/l");
    tOut->Branch("TimeStamp_ns_1", &TimeStamp_ns_1, "TimeStamp_ns_1/l");
    tOut->Branch("BeamGate_ns_0", &BeamGate_ns_0, "BeamGate_ns_0/l");
    tOut->Branch("BeamGate_ns_1", &BeamGate_ns_1, "BeamGate_ns_1/l");
    tOut->Branch("TimeStamp_ps_0", &TimeStamp_ps_0, "TimeStamp_ps_0/l");
    tOut->Branch("TimeStamp_ps_1", &TimeStamp_ps_1, "TimeStamp_ps_1/l");
    tOut->Branch("BeamGate_ps_0", &BeamGate_ps_0, "BeamGate_ps_0/l");
    tOut->Branch("BeamGate_ps_1", &BeamGate_ps_1, "BeamGate_ps_1/l");
    tOut->Branch("EventDeviation_ns_0", &EventDeviation_ns_0, "EventDeviation_ns_0/l");
    tOut->Branch("EventDeviation_ns_1", &EventDeviation_ns_1, "EventDeviation_ns_1/l");
    tOut->Branch("CTCTriggerIndex_0", &CTCTriggerIndex_0, "CTCTriggerIndex_0/l");
    tOut->Branch("CTCTriggerIndex_1", &CTCTriggerIndex_1, "CTCTriggerIndex_1/l");
    tOut->Branch("CTCTriggerTimeStamp_ns_0", &CTCTriggerTimeStamp_ns_0, "CTCTriggerTimeStamp_ns_0/l");
    tOut->Branch("CTCTriggerTimeStamp_ns_1", &CTCTriggerTimeStamp_ns_1, "CTCTriggerTimeStamp_ns_1/l");
    tOut->Branch("BGMinusTrigger_ns_0", &BGMinusTrigger_ns_0, "BGMinusTrigger_ns_0/L");
    tOut->Branch("BGMinusTrigger_ns_1", &BGMinusTrigger_ns_1, "BGMinusTrigger_ns_1/L");
    tOut->Branch("BGCorrection_tick_0", &BGCorrection_tick_0, "BGCorrection_tick_0/l");
    tOut->Branch("BGCorrection_tick_1", &BGCorrection_tick_1, "BGCorrection_tick_1/l");
    tOut->Branch("TSCorrection_tick_0", &TSCorrection_tick_0, "TSCorrection_tick_0/l");
    tOut->Branch("TSCorrection_tick_1", &TSCorrection_tick_1, "TSCorrection_tick_1/l");
    tOut->Branch("LAPPD_PPS_interval_ticks_0", &LAPPD_PPS_interval_ticks_0, "LAPPD_PPS_interval_ticks_0/l");
    tOut->Branch("LAPPD_PPS_interval_ticks_1", &LAPPD_PPS_interval_ticks_1, "LAPPD_PPS_interval_ticks_1/l");
    tOut->Branch("BG_PPSBefore_tick_0", &BG_PPSBefore_tick_0, "BG_PPSBefore_tick_0/l");
    tOut->Branch("BG_PPSBefore_tick_1", &BG_PPSBefore_tick_1, "BG_PPSBefore_tick_1/l");
    tOut->Branch("BG_PPSAfter_tick_0", &BG_PPSAfter_tick_0, "BG_PPSAfter_tick_0/l");
    tOut->Branch("BG_PPSAfter_tick_1", &BG_PPSAfter_tick_1, "BG_PPSAfter_tick_1/l");
    tOut->Branch("BG_PPSDiff_tick_0", &BG_PPSDiff_tick_0, "BG_PPSDiff_tick_0/l");
    tOut->Branch("BG_PPSDiff_tick_1", &BG_PPSDiff_tick_1, "BG_PPSDiff_tick_1/l");
    tOut->Branch("BG_PPSMissing_tick_0", &BG_PPSMissing_tick_0, "BG_PPSMissing_tick_0/l");
    tOut->Branch("BG_PPSMissing_tick_1", &BG_PPSMissing_tick_1, "BG_PPSMissing_tick_1/l");
    tOut->Branch("TS_PPSBefore_tick_0", &TS_PPSBefore_tick_0, "TS_PPSBefore_tick_0/l");
    tOut->Branch("TS_PPSBefore_tick_1", &TS_PPSBefore_tick_1, "TS_PPSBefore_tick_1/l");
    tOut->Branch("TS_PPSAfter_tick_0", &TS_PPSAfter_tick_0, "TS_PPSAfter_tick_0/l");
    tOut->Branch("TS_PPSAfter_tick_1", &TS_PPSAfter_tick_1, "TS_PPSAfter_tick_1/l");
    tOut->Branch("TS_PPSDiff_tick_0", &TS_PPSDiff_tick_0, "TS_PPSDiff_tick_0/l");
    tOut->Branch("TS_PPSDiff_tick_1", &TS_PPSDiff_tick_1, "TS_PPSDiff_tick_1/l");
    tOut->Branch("TS_PPSMissing_tick_0", &TS_PPSMissing_tick_0, "TS_PPSMissing_tick_0/l");
    tOut->Branch("TS_PPSMissing_tick_1", &TS_PPSMissing_tick_1, "TS_PPSMissing_tick_1/l");
    tOut->Branch("TS_driftCorrection_ns_0", &TS_driftCorrection_ns_0, "TS_driftCorrection_ns_0/l");
    tOut->Branch("TS_driftCorrection_ns_1", &TS_driftCorrection_ns_1, "TS_driftCorrection_ns_1/l");
    tOut->Branch("BG_driftCorrection_ns_0", &BG_driftCorrection_ns_0, "BG_driftCorrection_ns_0/l");
    tOut->Branch("BG_driftCorrection_ns_1", &BG_driftCorrection_ns_1, "BG_driftCorrection_ns_1/l");
    
    std::ofstream outputEvents("outputEvents.txt");
    for (auto it = ResultMap.begin(); it != ResultMap.end(); it++)
    {
        string key = it->first;
        vector<vector<ULong64_t>> Result = it->second;

        if (Result.size() == 0)
            continue;

        runNumber_out = std::stoi(key.substr(0, key.find("_")));
        subRunNumber_out = std::stoi(key.substr(key.find("_") + 1, key.find("_", key.find("_") + 1) - key.find("_") - 1));
        partFileNumber_out = std::stoi(key.substr(key.find("_", key.find("_") + 1) + 1, 
		key.find("_", key.find("_", key.find("_") + 1) + 1) - key.find("_", key.find("_") + 1) - 1));
        LAPPD_ID_out = std::stoi(key.substr(key.find("_", key.find("_", key.find("_") + 1) + 1) + 1, 
		key.size() - key.find("_", key.find("_", key.find("_") + 1) + 1) - 1));
        
        if (Result.size() < 48) {
            std::cerr << "Unexpected Result size = " << Result.size()
                  << " for key " << key << std::endl;
            continue;
        }

        final_offset_ns_0 = Result[0][0];
        final_offset_ns_1 = Result[24][0];
        final_offset_ps_negative_0 = Result[0][1];
        final_offset_ps_negative_1 = Result[24][1];
        gotOrphanCount_0 = Result[0][2];
        gotOrphanCount_1 = Result[24][2];
        gotMin_mean_dev_noOrphan_0 = Result[0][3];
        gotMin_mean_dev_noOrphan_1 = Result[24][3];
        increament_times_0 = Result[0][4];
        increament_times_1 = Result[24][4];
        min_mean_dev_0 = Result[0][5];
        min_mean_dev_1 = Result[24][5];
        EventNumberInThisPartFile = Result[0][8];
        drift0 = Result[0][9];
        drift1 = Result[24][9];

        //                                        0         1             2            3             4            5             6           
	//      7           8                9                10                       11                        12                          
	//       13                     14            15           16          17          18            19           20          21
        // vector<vector<ULong64_t>> Result = {FitInfo, TimeStampRaw, BeamGateRaw, TimeStamp_ns, BeamGate_ns, TimeStamp_ps, BeamGate_ps, 
	// EventIndex, EventDeviation_ns, CTCTriggerIndex, CTCTriggerTimeStamp_ns, BeamGate_correction_tick, TimeStamp_correction_tick, 
	// LAPPD_PPS_interval_ticks, BG_PPSBefore, BG_PPSAfter, BG_PPSDiff, BG_PPSMiss, TS_PPSBefore, TS_PPSAfter, TS_PPSDiff, TS_PPSMiss};

        // any Result[x] , if x > 13, x = x + 8
        for (int j = 0; j < Result[1].size(); j++)
        {
	    // 325250 ns is the software delay between CTC UBT and LAPPD Beamgate
            long long BGTdiff_0 = Result[4][j] - Result[10][j] - 325250; 
            long long BGTdiff_1 = Result[28][j] - Result[34][j] - 325250;

            // std::cout << "BGTDiff: " << BGTdiff << std::endl;
            // std::cout << "Saving BeamGate_ns = " << Result[4][j] << ", CTCTriggerTimeStamp_ns = " << Result[10][j]
		// << ", with BG-T-325250 = " << BGTdiff << ", at partFileNumber " << partFileNumber_out << ", EventIndex = "
		// << Result[7][j] << ", j = " << j << std::endl;
            
            outputEvents << fixed << Result[3][j] << " " << Result[5][j] << " " << Result[4][j] << " " 
		<< Result[6][j] << " " << Result[10][j] << " " << BGTdiff_0 << " " << partFileNumber_out 
		<< " " << Result[7][j] << " " << Result[11][j] << " " << Result[12][j] << " " << Result[13][j] 
		<< " " << Result[22][j] << " " << Result[23][j] << std::endl;
            
            EventIndex = Result[7][j];

            // --- Saving ACDC0 offsets ---
            TimeStampRaw_0              = Result[1][j];
            BeamGateRaw_0               = Result[2][j];
            TimeStamp_ns_0              = Result[3][j];
            BeamGate_ns_0               = Result[4][j];
            TimeStamp_ps_0              = Result[5][j];
            BeamGate_ps_0               = Result[6][j];
            CTCTriggerIndex_0           = Result[9][j];
            CTCTriggerTimeStamp_ns_0    = Result[10][j];
            EventDeviation_ns_0         = Result[8][j];
            BGMinusTrigger_ns_0         = BGTdiff_0;
            BGCorrection_tick_0         = Result[11][j];
            TSCorrection_tick_0         = Result[12][j];
            LAPPD_PPS_interval_ticks_0  = Result[13][j];
            BG_PPSBefore_tick_0         = Result[14][j];
            BG_PPSAfter_tick_0          = Result[15][j];
            BG_PPSDiff_tick_0           = Result[16][j];
            BG_PPSMissing_tick_0        = Result[17][j];
            TS_PPSBefore_tick_0         = Result[18][j];
            TS_PPSAfter_tick_0          = Result[19][j];
            TS_PPSDiff_tick_0           = Result[20][j];
            TS_PPSMissing_tick_0        = Result[21][j];
            TS_driftCorrection_ns_0     = Result[22][j];
            BG_driftCorrection_ns_0     = Result[23][j];

            // --- Saving ACDC1 offsets --- 
            TimeStampRaw_1              = Result[25][j];
            BeamGateRaw_1               = Result[26][j];
            TimeStamp_ns_1              = Result[27][j];
            BeamGate_ns_1               = Result[28][j];
            TimeStamp_ps_1              = Result[29][j];
            BeamGate_ps_1               = Result[30][j];
            CTCTriggerIndex_1           = Result[33][j];
            CTCTriggerTimeStamp_ns_1    = Result[34][j];
            EventDeviation_ns_1         = Result[32][j];
            BGMinusTrigger_ns_1         = BGTdiff_1;
            BGCorrection_tick_1         = Result[35][j];
            TSCorrection_tick_1         = Result[36][j];
            LAPPD_PPS_interval_ticks_1  = Result[37][j];
            BG_PPSBefore_tick_1         = Result[38][j];
            BG_PPSAfter_tick_1          = Result[39][j];
            BG_PPSDiff_tick_1           = Result[40][j];
            BG_PPSMissing_tick_1        = Result[41][j];
            TS_PPSBefore_tick_1         = Result[42][j];
            TS_PPSAfter_tick_1          = Result[43][j];
            TS_PPSDiff_tick_1           = Result[44][j];
            TS_PPSMissing_tick_1        = Result[45][j];
            TS_driftCorrection_ns_1     = Result[46][j];
            BG_driftCorrection_ns_1     = Result[47][j];

            tOut->Fill();
        }
        
        outputOffset << runNumber_out << "\t" << subRunNumber_out << "\t" << partFileNumber_out << "\t" << 0 
		<< "\t" << LAPPD_ID_out << "\t" << final_offset_ns_0 << "\t" << final_offset_ns_1 << "\t" 
		<< final_offset_ps_negative_0 << "\t" << final_offset_ps_negative_1 << "\t" << gotOrphanCount_0 
		<< "\t" << gotOrphanCount_1 << "\t" << EventNumberInThisPartFile << "\t" << gotMin_mean_dev_noOrphan_0 
		<< "\t" << gotMin_mean_dev_noOrphan_1 << "\t" << increament_times_0 << "\t" << increament_times_1 
		<< "\t" << min_mean_dev_0 << "\t" << min_mean_dev_1 << std::endl;
    }

    outputOffset.close();
    fOut->cd();
    tOut->Write();
    fOut->Close();
    cout << "Result saved." << endl;
}
