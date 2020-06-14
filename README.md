# PXD013277_E-coli_spike-ins_MS2-TMT

E. coli spiked into a human background, 10-plex TMT labeled, MS2 scans, Q Exactive

## Overview

This is data set D3 (E. coli Spike-In TMT Data Set) from [PXD013277](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD013277) PRIDE archive from [this publication](https://www.mcponline.org/content/19/6/1047):

> Zhu, Y., Orre, L.M., Tran, Y.Z., Mermelekas, G., Johansson, H.J., Malyutina, A., Anders, S. and Lehtiö, J., 2020. DEqMS: a method for accurate variance estimation in differential protein expression analysis. Molecular & Cellular Proteomics, 19(6), pp.1047-1057.

This is a common type of test mixture where one complex proteome (E. coli lysate) is spiked into a complex background (MCF7 human cell lysates) in a series of known relative amounts. The E. coli lysate was spiked in at 3 different relative levels as indicated in the table below:

Channel|E. coli (ug)|MCF7 (ug)|Total (ug)|Human Fraction
-------|------------|---------|----------|--------------
126C|7.5|70|77.5|0.90
127N|7.5|70|77.5|0.90
127C|7.5|70|77.5|0.90
128N|15|70|85|0.82
128C|15|70|85|0.82
129N|15|70|85|0.82
129C|15|70|85|0.82
130N|45|70|115|0.61
130C|45|70|115|0.61
131N|45|70|115|0.61

The mixes were labeled with [tandem mass tag](https://pubs.acs.org/doi/abs/10.1021/ac0262560) (TMT) 10-plex reagents and analyzed in a large scale fractionation using a Q Exactive HF instrument. Briefly, samples were reduced and alkylated, digested with LysC and trpysin in an SP3 bead protocol. The first dimension separation was a peptide IEF-IPG method that resulted in 72 fractions. Each fraction was run on a standard reverse phase column into a Q Exactive HF instrument. Instrument parameters looked typical for an MS2 TMT experiment on this platform.

The human samples and the E. coli samples are replicates of one cell lysate pellet (I think), and the digestions and labeling reactions are all independent. We will have some sample processing variation, but not as much as true independent biological replicates. One sentence in the Methods section (bottom of first column of page 1049) caught my eye:

> "The statistical analysis was focused on the comparison between the 3 replicates with 7.5 microgram and the 4 replicates with 15 microgram of spiked in E. coli protein extract (the group with 45 microgram spike-in had substantial global impact on protein ratios after normalization)."

This looked like a pretty classic two proteome test mixture experiment. Human and E. coli have low tryptic peptide sequence homology, so issues of shared peptides are minimal. Why would the 45 microgram channels not be usable?

My [PAW pipeline](https://github.com/pwilmart/PAW_pipeline) is the best tool for understanding TMT data sets. Protein expression is kept in the natural reporter ion peak height intensity scale. Sensitive and accurate peptide spectrum match (PSM) identifications are obtained using Comet and the [target/decoy method](https://www.nature.com/articles/nmeth1019?report=reader). Basic and extended parsimony protein grouping is used to maximize usable quantitative information. Reporter ions from individual PSMs are filtered to remove very low signals, summed up into protein total intensities, and any final zero intensities replaced with a low signal (value of 150). The protein total intensities are good proxies for relative protein abundances.  

The data was downloaded from PRIDE (about 29 GB) and converted to compressed text files using [MSConvert](http://proteowizard.sourceforge.net/). The re-analysis used the [PAW pipeline](https://link.springer.com/article/10.1007/s12177-009-9042-6). A combined human and E. coli canonical reference proteome FASTA file was created using [fasta_utilities tools](https://github.com/pwilmart/fasta_utilities). A wide tolerance [Comet search](http://comet-ms.sourceforge.net/) (1.25 Da for parent ions) was used with recommended settings for high-resolution MS2 fragment ions. PSMs were filtered to 1% FDR. Protein inference used basic and extended parsimony analysis. TMT analysis used reporter ions from all unique (w.r.t. the final protein list) peptides summed into protein total intensities. Normalizations and statistical testing was done in [edgeR](https://academic.oup.com/bioinformatics/article/26/1/139/182458) using [Jupyter notebooks](https://jupyter.org/).

---

## PAW Pipeline FDR Analysis

![Deltamass 2+](images/deltamass_2plus.png)

2+ peptide deltamasses. We have relatively minor amounts of deamidation (peak at 0.984 Da). Not sure if the -1 Da peak has any significance. The MS2 spectra are high resolution and Comet is configured with a small 0.02 Da fragment ion tolerance. Comet results look a little different for high res MS2 datasets.

---

![Deltamass 3+](images/deltamass_3plus.png)

3+ peptide deltamasses look very similar to 2+ peptides. We probably have about twice as many 2+ peptides as 3+ peptides.

---

![Deltamass 4+](images/deltamass_4plus.png)

These are the 4+ peptide deltamass distributions. Interestingly, the -1 Da peak is gone. We have more first isotopic peak triggers (peak at 1.003 Da) for higher charge state peptides. We can also see that the deltamass peak widths increase with increasing charge.

---

![Scores 2+ 0-Da nomods](images/scores_0Da_2plus_nomod.png)

Unmodified 2+ peptides for the 0-Da deltamass window. The target matches are the blue distributions and the decoy matches are the red distributions. The dotted line is the 1% FDR cutoffs.

---

![Scores 2+ 0-Da oxmet](images/scores_0Da_2plus_oxmet.png)

Oxidized Met 2+ peptides for the 0-Da deltamass window.

---

![Scores 3+ 0-Da nomods](images/scores_0Da_3plus_nomod.png)

Unmodified 3+ peptides for the 0-Da deltamass window.

---

![Scores 3+ 0-Da oxmet](images/scores_0Da_3plus_oxmet.png)

Oxidized Met 3+ peptides for the 0-Da deltamass window.

---

![Scores 4+ 0-Da nomods](images/scores_0Da_4plus_nomod.png)

Unmodified 4+ peptides for the 0-Da deltamass window.

---

![Scores 4+ 0-Da oxmet](images/scores_0Da_4plus_oxmet.png)

Oxidized Met 4+ peptides for the 0-Da deltamass window.

---

![Scores 2+ 1-Da nomods](images/scores_1Da_2plus_nomod.png)

Unmodified 2+ peptides for the 1-Da deltamass window.

---

![Scores 2+ 1-Da oxmet](images/scores_1Da_2plus_oxmet.png)

Oxidized Met 2+ peptides for the 1-Da deltamass window.

---


![Scores 3+ 1-Da nomods](images/scores_1Da_3plus_nomod.png)

Unmodified 3+ peptides for the 1-Da deltamass window.

---

![Scores 3+ 1-Da oxmet](images/scores_1Da_3plus_oxmet.png)

Oxidized Met 3+ peptides for the 1-Da deltamass window.

---

![Scores 4+ 1-Da nomods](images/scores_1Da_4plus_nomod.png)

Unmodified 4+ peptides for the 1-Da deltamass window.

--

![Scores 4+ 1-Da oxmet](images/scores_1Da_4plus_oxmet.png)

Oxidized Met 4+ peptides for the 1-Da deltamass window.

---

![Scores 2+ outside nomods](images/scores_noDa_2plus_nomod.png)

Unmodified 2+ peptides with deltamasses outside of 0-Da and 1-Da windows.

---

![Scores 2+ outside oxmet](images/scores_noDa_2plus_oxmet.png)

Oxidized Met 2+ peptides with deltamasses outside of 0-Da and 1-Da windows.

---

![Scores 3+ outside nomods](images/scores_noDa_3plus_nomod.png)

Unmodified 3+ peptides with deltamasses outside of 0-Da and 1-Da windows.

---

![Scores 3+ outside oxmet](images/scores_noDa_3plus_oxmet.png)

Oxidized Met 3+ peptides with deltamasses outside of 0-Da and 1-Da windows.

---

![Scores 4+ outside nomods](images/scores_noDa_4plus_nomod.png)

Unmodified 4+ peptides with deltamasses outside of 0-Da and 1-Da windows.

---

![Scores 4+ outside oxmet](images/scores_noDa_4plus_oxmet.png)

Oxidized Met 4+ peptides with deltamasses outside of 0-Da and 1-Da windows.

---

The RAW data files had 1,961,506 MS2 scans. At a 1% FDR, there were 303,201 scans that passed the score cutoffs (15.5% overall ID rate). There were 9,650 protein IDs (using a two peptides per protein criterion) excluding common contaminants, decoys, and proteins without any usable reporter ion signals.

---

## High-level Quantitative Overview

The total reporter ion intensity for each channel is a good proxy for the total amount of protein in that sample. The protein totals in the PAW pipeline exclude any shared peptides. However, homologous protein grouping greatly reduces the numbers of shared peptides in the results summaries and maximizes the number of usable reporter ion signals for quantification.

> The PAW pipeline uses reporter ion peak heights for its quantitative values. These are well correlated with total protein MS1 feature intensities in MaxQuant. Our group also studies human and animal eye lenses. It has been well established that a small number of lens specific crystallin proteins make up over 90% of the proteins in lens. MS1 feature intensities correctly estimate that the lens crystallins are more than 90% of the total intensity. Spectral counting estimates that about 30% of the counts are from lens crystallins. It is well known that spectral counting underestimates high abundance proteins and overestimates low abundance proteins. Summed reporter ions in TMT experiments also estimates that over 90% of the signal in lens comes from crystallins. Reporter ion signal-to-noise ratios (the default setting in Proteome Discoverer) **incorrectly** estimate that about 60% of the lens proteins are crystallins. Signal-to-noise ratios, a compressed and unitless quantity, are not valid abundance measures.

We can take the column totals (after excluding common contaminants and decoy matches) for each channel to see what the protein totals across the channels look like.  

Channel|Spike-in|All|Human|E. coli
-------|--------|---|-----|-------
126C|7.5|103,291,224,641|86,089,119,909|17,202,104,732
127N|7.5|94,677,244,028|79,728,913,187|14,948,330,840
127C|7.5|100,610,569,631|83,563,004,123|17,047,565,507
128N|15|120,468,775,251|91,414,111,974|29,054,663,276
128C|15|115,258,220,367|90,022,292,892|25,235,927,475
129N|15|126,712,117,536|97,032,335,772|29,679,781,764
129C|15|104,602,264,635|81,434,898,939|23,167,365,696
130N|45|110,693,888,205|69,036,827,572|41,657,060,633
130C|45|122,603,435,090|75,677,201,918|46,926,233,171
131N|45|114,637,776,618|70,082,564,633|44,555,211,985

![All Totals](images/All_totals.png)

We see that there is some variation in the protein totals across the channels. Labeling protocols usually start with the same total amount of proteins. We usually see variation like this in TMT data sets. It always seems larger than it should be. The x-axis indicates the micrograms of E. coli lysate added to the 70 micrograms of human background. The total amount of human + E. coli protein was the same in each channel.

![E. coli Totals](images/E-coli_totals.png)

From the protein accessions, we can tell which proteins are human and which are E. coli. We can get the column totals for each species. These are the E. coli intensity totals. We see increased intensity totals as the amount of E. coli lysate is increased. However, the relative changes in intensities are about half the expected values (the 45s are about 3 times the 7.5s instead of the expected 6 times).

![Human Totals](images/Human_totals.png)

The human proteins are not too different except for the last three where the larger amount of E. coli pushes the human background down a little.

We can average the intensities for each mixture and see how the total amounts of protein compare.

Spike-In Level|Average Human Total|Average E. coli Total|Expected Fraction|Actual Fraction
--------------|-------------------|---------------------|-----------------|---------------
7.5 out of 70|83.1 Billion|16.4 Billion|11%|20%
15 out of 70|90.0 Billion|26.8 Billion|21%|30%
45 out of 70|71.6 Billion|44.4 Billion|64%|62%

E. Coli Levels:

Levels|Expected Ratio|Observed Ratio
------|--------------|--------------
7.5 to 15|2|1.6
7.5 to 45|6|2.7
15 to 45|3|1.7

It has been known for several years that reporter ion accuracy in MS2-based isobaric labeling seems rather poor. Precision, however, seems quite good. The situation is very different in the SPS MS3 methods available on newer Tribrid Orbitraps. We had accurate reporter ion intensities in the dilution series experiment (see: https://github.com/pwilmart/Dilution_series) using SPS MS3. Both intensity column totals and individual protein expression trend lines were consistent with the known dilution factors.

---

## Individual Protein Measurements
