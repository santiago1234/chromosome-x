Link autosomes

```bash
for i in {1..22}; do
  ln -s /home/smedina/results/chromosome-x/250422-raw-vcfs/tgp_chr${i}.vcf.gz tgp_chr${i}.vcf.gz
done
```

We now link the processed chromosome X.

```bash
ln -s /home/smedina/results/chromosome-x/250423-processed_chrX/chrX.vcf.gz tgp_chrX.vcf.gz
```
