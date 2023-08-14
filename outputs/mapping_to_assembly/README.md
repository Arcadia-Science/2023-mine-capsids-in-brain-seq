Code to assemble viral genomes from mapped reads.
```
mkdir mapping_to_assembly
cd mapping_to_assembly
ln -s ../mapping/SRR1778915_nucleotides_pe.sorted.bam .
ln -s ../mapping/SRR1779200_nucleotides_pe.sorted.bam .
ln -s ../mapping/SRR1779200_nucleotides_pe.sorted.bam.bai .
ln -s ../mapping/SRR1778915_nucleotides_pe.sorted.bam.bai .
mamba install bedtools bbmap megahit=1.2.8
bamToFastq -i SRR1778915_nucleotides_pe.sorted.bam -fq SRR1778915_nucleotides_pe.fq
repair.sh in=SRR1778915_nucleotides_pe.fq out=SRR1778915_nucleotides_pe_R1.fq out2=SRR1778915_nucleotides_pe_R2.fq outs=SRR1778915_nucleotides_pe_orphan.fq repair=t overwrite=true
megahit -1 SRR1778915_nucleotides_pe_R1.fq -2 SRR1778915_nucleotides_pe_R2.fq -o SRR1778915
# mv to a better file name
mv SRR1778915/final.contigs.fa SRR1778915_final.contigs.fa
rm -rf SRR1778915
```

deleted these files since the assemblies didn't add anything
```
bamToFastq -i SRR1779200_nucleotides_pe.sorted.bam -fq SRR1779200_nucleotides_pe.fq
repair.sh in=SRR1779200_nucleotides_pe.fq out=SRR1779200_nucleotides_pe_R1.fq out2=SRR1779200_nucleotides_pe_R2.fq outs=SRR1779200_nucleotides_pe_orphan.fq repair=t overwrite=true
megahit -1 SRR1779200_nucleotides_pe_R1.fq -2 SRR1779200_nucleotides_pe_R2.fq -o SRR1779200


# try a combined assembly
cat SRR1779200_nucleotides_pe_R1.fq SRR1778915_nucleotides_pe_R1.fq > R1.fq
cat SRR1779200_nucleotides_pe_R2.fq SRR1778915_nucleotides_pe_R2.fq > R2.fq
megahit -1 R1.fq -2 R2.fq -o combined
```
