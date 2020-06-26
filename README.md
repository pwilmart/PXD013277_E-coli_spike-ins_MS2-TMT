# PXD013277_E-coli_spike-ins_MS2-TMT

E. coli spiked into a human background, 10-plex TMT labeled, MS2 scans, Q Exactive

- [Data Overview](#overview)
- [PAW Pipeline](#pipeline)
- [High-level Quantitities](#high-level)
- [Human Proteins](#human)
- [Explanation of distortion](#theory)
- [E. coli Proteins](#E.coli)
- [Summary](#summary)

## <a name="overview"></a>Overview

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

## <a name="pipeline"></a>PAW Pipeline FDR Analysis

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

## <a name="high-level"></a>High-level Quantitative Overview

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

## <a name="human"></a>Human Protein Measurements

The established way that these spike-in experiments are analyzed is to focus on the expression levels of the E. coli proteins and completely ignore the background human proteome. That is insufficient. We always have to normalize the data. We also usually do some differential expression testing. The behavior of the background of unchanging proteins is critical for both normalization and statistical testing. In this study design, the human proteome is the constant background.

The notebooks I used for this data are derived from the typical notebooks I use for TMT analysis. I like using the Bioconductor package edgeR. That package has a nice normalization routine called the [trimmed mean of M-values](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25) (TMM). This routine removes some of the highest and lowest abundance data points. It computes the M values (fold changes) and trims those to remove large and small fold changes. This should define a core set of proteins that are the unchanged proteins. The algorithm calculates the single multiplicative factors that make those proteins as similar as possible.

Three normalization methods were used in these notebooks:

- [TMM normalization on all 9,650 proteins](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons.html)
- [Manual normalization to equalize human protein total intensity](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons_no-norm.html)
- [Filter dataset for human proteins only, then run TMM](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons_human.html)

The differential expression testing used edgeR and DE candidates were proteins with a expression FDR (Benjamini-Hochberg corrected statistical model p-values) of 10% or less (FDR < 0.10). In the `TMM ALL` and the `Manual All` methods, the testing was done on all 9,650 proteins. The `TMM Human` used just the 7,559 human proteins. The E. coli proteins are only over-expressed relative to the human proteins (that is what we mean by spiked in). That is sort of a situation that the TMM algorithm was designed for. The trimming might exclude most of the E. coli proteins and then the TMM factors would be making the human background more similar.

If there are too many DE proteins, the algorithm might not work well. How many E. coli proteins can we have and still have TMM work okay? The sums of all of the human protein reporter ions can be computed for each channel. Factors to match each channel to have the same total for the human proteins can be computed manually. We can do the edgeR analysis without using the TMM function (`Manual All` method).

If there is the chance that the presence of the E. coli proteins could throw off the the TMM algorithm, we can filter out the E. coli proteins and just work with the human background proteome in each channel. We can then run the TMM function and expect that the human background will be normalized as well as possible.

The table below shows how many DE candidates the statistical testing falsely identifies. The human background is the same 70 microgram of MCF7 cell lysate for each channel. The ground truth is that all human proteins should be at the same relative intensities after we have properly normalized the data. The last line should be 7,559 across the row. We see a large difference as a function of the amount of E. coli that was spiked in. We also have a dependence on the normalization method.  

<br />  

Method|TMM All|Manual All|TMM Human|TMM All|Manual All|TMM Human|TMM All|Manual All|TMM Human
------|-------|----------|---------|-------|----------|---------|-------|----------|---------   
Conditions|7.5 vs 15|7.5 vs 15|7.5 vs 15|15 vs 45|15 vs 45|15 vs 45|7.5 vs 45|7.5 vs 45|7.5 vs 45
E. coli difference|1:2|1:2|1:2|1:3|1:3|1:3|1:6|1:6|1:6
All Human Proteins|7,559|7,559|7,559|7,559|7,559|7,559|7,559|7,559|7,559
DE Up|151|256|81|525|1,275|624|675|1,766|956
DE Down|228|47|6|1,644|161|245|2,549|359|688
non-DE|7,180|7,256|7,472|5,390|6,123|6,690|4,335|5,434|5,915

<br />

We can also look at the data from the last three rows in the table as a bar plot. Clearly, the E. coli proteins are somehow affecting both the statistical testing and the normalization. We have the lowest levels of false human DE candidates when we exclude the E. coli proteins and then run TMM. We still have an even larger effect associated with the amount of E. coli that was spiked in. I have Jupyter notebooks where all of this testing was done.

* [All proteins with TMM](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons.html)
* [All proteins without TMM](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons_no-norm.html)
* [Human proteins with TMM](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons_human.html)

> Remember that the human proteins are a constant background in all 10 channels. They should be the same and no human proteins should be differentially expressed.

<br />


![Human candidates](images/human_candidates.png)

<br />

We can see how well the TMM normalization worked for the human proteins by looking at some distribution boxplots. We have very good horizontal alignment of the medians (the notches) and the interquartile ranges (the actual boxes). This global view of the normalized data looks ideal. It looks like the human proteins are the same in all channels.

![Boxplots](images/boxplots.png)

We can also check the clustering of the normalized data. We see that the human proteins strongly cluster by the amount of E. coli that was spiked in. This supports the idea that the E. coli proteins themselves are altering the reporter ion intensities of the human proteins. We do not actually have the human proteome background being the same in each channel. The human expression is distorted by the E. coli proteins.

![Clustering](images/MDS_plot.png)

This is the MA plot for just the human proteins (with TMM norm) for the comparison of the 7.5 microgram and 15 microgram spike-in channels. We mostly have non-DE candidates (as expected) for the human proteins. We have a small number of statistically significant candidates. The number is not large, but the candidates are biased towards over-expression in the samples where more E. coli was present.

![Not so much E. coli](images/human_7.5_15.png)

<br />

Here is the MA plot for the 7.5 microgram and 45 microgram spike-in channels. Something has gone really wrong! The blue points all look like false positives. The red points are more complicated. Some seem to be a mirror image of the blue false positives. Others seem to be truly over-expressed. This is a real results train wreck.

![Lots of E. coli](images/human_7.5_45.png)

---

## <a name="theory"></a>What is happening?

We have 1,644 DE candidates from the unchanged human background proteins! That is almost 22% of the proteins. We have done a high quality normalization. We have a test statistic using a moderated test statistic. We have done multiple testing corrections.

Here is my **theory**: I think there is a low level background of previously sampled peptides that is coming off of the column. Kind of a "ghost of Christmas past" situation (the ghost of peptides past...). A precursor might have 3 to 4 very narrow isotopic peaks. The isolation quadrupole might be 2 Da on older instruments and maybe 0.7 Da on newer platforms. The isolation window is 2000 milliDa to 700 milliDa wide compared to a few peaks that may be 10 milliDa wide at the base. Low levels of background signal can be integrated over a very wide range compared to the widths of the precursor peaks.

Distortions in isobaric labeling reagent tags acquired in MS2 scans has been known for years and documented in many published papers. When high resolution instruments became the norm, it was easier to look for any other co-isolated precursors based on what fraction of the isolated peaks were mapped to the precursor of interest. These precursor purity filters looked great in theory; however, they did not fix the reporter ion distortions. Along the same lines, larger-scale peptide separations should reduce the chance of co-eluting peptides compared to less separation. That also failed to improve the ratio compression problem.

> Keep in mind that we have a very large separation scheme for this data. It is something like 72 fractions. The PAW pipeline does not have any way to check for precursor purity (it does not use any information from survey scans).

I think the only explanation left is the non-specific background effect. I wrote a blog about how that explains ratio compression (https://pwilmart.github.io/blog/2020/01/05/TMT-ratio-distortions). This spike-in experiment actually provides an excellent test of my theory. We have a similar human background proteome for all channels. If we get some other human proteins in the non-specific noise, they will not change the expression pattern of the human protein of interest. However, any E. coli proteins that sneak through into the non-specific background do not have an expression pattern like the human background. That is the whole point, right? We have a small amount of E. coli in the first 3 channels and 6 times as much E. coli in the last three channels. Any contribution to the reporter ions from any E. coli proteins will boost the signal for the human protein in the last three channels. For the 45 microgram E. coli spike-in, the E. coli is about 61% of the total protein. We have over 2,000 E. coli proteins out of the 9,500 total proteins (about half of the proteome in both cases).

It is likely we have some boost in the last three channels for many of the human proteins due to some "ghosts of E. coli peptides past". How will this manifest itself? What sneaks through in the background will probably depend on the m/z and the retention time. Not all human proteins in the 45 microgram E. coli spike-in will be boosted. Some will and some won't. Our normalization methods typically assume most proteins do not change and try to make that non-differentially expressed background more similar between samples. These proteome spike-in experiments are not very compatible with that logic because the spike-in can only go in one direction. The differences are just how much all of the spike-in proteins are over expressed.

If we have some bleed through of E. coli peptides and have some human peptides get some partial reporter ion intensity boosts, then the total signal of the human proteins (our background here) for the 45 microgram E. coli spike-in will be larger than it really should have been. This will alter the single global scaling factor. The human proteins in the 45 microgram E. coli spike-in have some boosted and some not boosted. The global normalization factor will end up splitting the difference. The boosted proteins will still be larger than they should have been. The true background level would not have been estimated correctly and the non-boosted proteins will not be fully brought up to the proper background levels. Some of the non-changed human proteins in the 45 microgram E. coli spike-in will end up with lower intensities that they should have been after this split-the-difference scaling.

I think that effect is pretty obvious in the MA plot above. All of the blue proteins (reduced expression) should be black and a mirror group of red proteins should also be unchanged (black). The more heavily over-expressed (red) proteins are the human proteins that had more E. coli boost.

We have a **double whammy effect** going on. Individual protein expression levels can be altered by the ghost background. In most experiments, we do not have anything like we have here with the high spike-in of E. coli. The background proteins are more likely to be "flat" (not differentially expressed), and this biases expression ratios to the mean expression (larger fold-changes become smaller). The other effect is the collective distortion of the unchanged background levels of the non-differentially expressed proteins. That can mess with the normalization methods and create issues for the statistical testing. Most statistical testing assumes that data are normalized correctly such that unchanged proteins would have similar means between conditions. We can get large numbers of false positive and false negative test results when we have unintentionally poorly normalized data.

---

## <a name="E.coli"></a>E. coli proteins

Most E. coli proteins were correctly determined to be differentially expressed in the three comparisons. We did have a few in the 2-fold comparison (about 10% of the proteins) that were false negatives. We do not get the expected ratios. The individual proteins generally matched the number we had above for the ratios of the grand totals of E. coli proteins.

<br />

Method|TMM All|Manual All|TMM All|Manual All|TMM All|Manual All
------|-------|----------|-------|----------|-------|----------
Conditions|7.5 vs 15|7.5 vs 15|15 vs 45|15 vs 45|7.5 vs 45|7.5 vs 45
E. Coli Proteins|2091|2091|2091|2091|2091|2091
DE Up|1883|1939|2062|2075|2072|2082
DE Down|0|0|0|0|1|0
non-DE|208|152|29|16|18|9

---

## <a name="summary"></a>Summary

It is clear that reporter ions from MS2 scans are not reliable measurements. In addition to the well established ratio compression effects that affect accuracy, we can have other data distortions that are impossible to correct. This experiment has a large number of differential proteins (all of the E. coli) in a specific expression pattern. This affects both normalization and statistical testing (partly because of its effect on normalization). Are most real experiments similar to this test mixture? Probably not. I do think that real world samples might still have potential for a lot of false positives from the unchanged background proteins. That makes for noisier DE candidate lists. That is not always a big problem if sufficient biological filtering can be done. Many time prior knowledge and biological patterns can be sued to filter out the noise.

I do not like the characteristics of any of the MS2-based TMT data I have looked at. I do think that SPS MS3 data is vastly superior. There are still dynamic range issues from isotopic carryover from adjacent channels that affect what types of quantitative measurements you want to use TMT to do. I do not think it is wise to attempt corrections for reagent purities. I do not think the measured values are all that accurate. N- and C-forms of the tags greatly complicate the corrections. Measurements at the PSM level are too noisy for the linear algebra corrections. Corrections of protein level summaries might be possible, but one has to take great care to do the aggregations, corrections, and normalizations in the proper order. I do not think the corrections are worth the effort.

I think TMT on platforms like the Q Exactives or the new Exploris line should not be done. I think the SPS MS3 method on Tribrids is the only viable quantitative TMT method. Honestly, I have nothing against the Exploris or Q Exactive instruments. They just are not capable of doing trustworthy TMT measurements.

---

#### Notebooks are here:

* [PXD013277_comparisons_human](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons_human.html) - Notebook looking at unchanged human background
* [PXD013277_comparisons_no-norm](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons_no-norm.html) - Notebook with manually matched human protein levels between spike-in channels
* [PXD013277_comparisons](https://pwilmart.github.io/PXD013277_E-coli_spike-ins_MS2-TMT/PXD013277_comparisons.html) - Notebook with generic edgeR workup (TMM normalization and exact testing)

---

Phil Wilmarth<br />June 16, 2020
