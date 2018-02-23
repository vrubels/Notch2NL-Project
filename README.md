# Notch2NL-Project
Repository for contributing scripts and explanatory figures around Notch2NL computational workflow


This walk through starts at the conclusion of the Gordion assembler pipeline, which constructs K assembled N2NL paratypes from 10X Chromium reads. In Chimp, we see a total of 19 assembled paratypes to start, which look like this against hg38.

![assembly](https://github.com/vrubels/Notch2NL-Project/blob/vrubels-edit-readme/Screen%20Shot%202018-02-23%20at%201.05.13%20PM.png)

Next, we step through known N2NL-associated Chimp BAC and other published sequences, to see if any of these assembled bins clearly support one or more known sequences. RNA fusion transcripts are included where available, as are annotations for any important genes we'd expect to see in this region. In general, we expect the assembler to have made 2 almost identical bins for each true locus, one each for the maternal and paternally derived versions of the region.

The sequence mappings I'll show below are:

* BACs
  * CH251-414E24
  * CH251-305M22
  * CH251-397N13
  * CH251-243D3
  * CH251-485H13
* Genomic Sequences
  * AACZ03011046.1
  * AADA01172734.1
  * AACZ03149796.1
  * AADA01170825.1
  * AACZ03149820.1
  * AADA01102044.1
 
 First for CH251-414E24, we see that c7 and c16 are the bins that most closely match the BAC reference.
 
 ![CH251-414E24](https://github.com/vrubels/Notch2NL-Project/blob/vrubels-edit-readme/Screen%20Shot%202018-02-23%20at%201.23.17%20PM.png)
 
 Next for CH251-305M22, we see that c0 and c18 are the bins that most closely match the BAC reference. 
 
 ![CH251-305M22](https://github.com/vrubels/Notch2NL-Project/blob/vrubels-edit-readme/Screen%20Shot%202018-02-23%20at%201.54.16%20PM.png)
 
 We see that CH251-397N13 has best evidence from the same c0 and c18 bins, so these BACs are grouped together.
 
 ![CH251-397N13](https://github.com/vrubels/Notch2NL-Project/blob/vrubels-edit-readme/Screen%20Shot%202018-02-23%20at%201.55.05%20PM.png)
 
 Finally, for CH251-243D3 and CH251-485H13, only c10 seems to be the best assembled bin. This is possibly due to a failure to separate similar maternal and paternal sequences into different bins, but no further bins are created upon reassembly.
 
 ![CH251-243D3](https://github.com/vrubels/Notch2NL-Project/blob/vrubels-edit-readme/Screen%20Shot%202018-02-23%20at%201.59.32%20PM.png)
 ![CH251-485H13](https://github.com/vrubels/Notch2NL-Project/blob/vrubels-edit-readme/Screen%20Shot%202018-02-23%20at%202.00.05%20PM.png)
 
 

The following table summarizes the results from this walk through, associating assembled bins with their names as they appear in our manuscript: https://doi.org/10.1101/221226

Name | Bin | Associated Genes | Genomic/BAC evidence
-----|-----|------------------|---------------------
delE2-only |	c13 | |
delE2-lrig |	c7, c16 | LRIG2 | CH251-414E24
delE2-sort	| c0, c18 | SORT1 | CH251-305M22/CH251-397N13
delE2-lrig+sort |	c10 | LRIG2, SORT1 | CH251-243D3/CH251-485H13
magi 3'	| c1, c4 | MAGI3 | AACZ03011046.1/AADA01172734.1
magi 5'	| c5, c14 | MAGI3 | AACZ03149796.1/AADA01170825.1
txnip	| c6, c8 | TXNIP | AACZ03149820.1
pde4dip	| c3, c15 | PDE4DIP | AADA01102044.1
notch2	| c2, c12, c17 | |
