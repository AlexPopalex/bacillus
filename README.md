This is the repository for the scripts used in the manuscript

### Evolution of alternative reproductive systems in Bacillus stick insects

by Guillaume Lavanchy*, Alexander Brandt*, Marc Bastardot, Zoé Dumas, Marjorie Labédan, William Toubiana, Patrick Tran Van, Andrea Luchetti, Valerio Scali, Barbara Mantovani, Tanja Schwander.



# mtDNA phylogeny
We added one reference individual per haplotype from this paper: https://doi.org/10.1006/mpev.2000.0850
We aligned the sequences with muscle 3.8.1551:
```
muscle -in All_filtered_with_refs.fa -phys > phylip.tmp
NSEQ=$( head -1 phylip.tmp | cut -f 1 -d ' ' )
NSITES=$( head -1 phylip.tmp | cut -f 2 -d ' ' )
echo $NSEQ $NSITES > All_filtered_with_refs_aligned.phylip
tail -n +2 phylip.tmp | tr -d '\n' | tr '@' '\n' | sed 's/ /@/1' | tr -d ' ' | tr '@' ' ' >> All_filtered_with_refs_aligned.phylip
```

The resulting alignment is available in this repository: `All_filtered_with_refs_aligned.phylip`. We then reconstructed the phylogeny with iqtree 2.0.6:
```
iqtree2 -T 1 -s Sicily_with_refs_aligned.phylip -m MFP -B 1000
```
The resulting tree was analysed in R along with the rest.
