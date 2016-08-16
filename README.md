# reco-ntuples
A modification off of the original Ntuplizer for the HGCAL reconstruction software studies.
This modification takes only the Rec Hits information needed for the deep-learning studies.

to be added on top of CMS-HGCAL:emilio_plus_sharing (the other two cms-merge-topics are temporary until https://github.com/CMS-HGCAL/cmssw/pull/2 and https://github.com/CMS-HGCAL/cmssw/pull/4 are merged):

```
cmsrel CMSSW_8_1_0_pre8
cd CMSSW_8_1_0_pre8/src
cmsenv
git cms-merge-topic CMS-HGCAL:emilio_plus_sharing
git cms-merge-topic lgray:topic_simclustering_pre7
git cms-merge-topic clelange:RecHitToolsImprovements
git clone https://github.com/kchang2/reco-ntuples RecoNtuples
scram b -j9
cd RecoNtuples/HGCalAnalysis/test
cmsRun twogamma_pt5_eta2_nosmear_calib.py
```

The input file needs to be step3 (i.e. RECO).
