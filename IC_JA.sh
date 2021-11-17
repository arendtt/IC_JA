# Os arquivos que você vai usar estão nessa pasta abaixo no Cerebro. São doze arquivos INPD_GenoBroad19_* . Cada um pesa 6.1Gb zipado. Já transferi o primeiro para a sua pasta

/genetica_1/Genotipagem_INPD_Broad2019/PGC_NIMH_Salum_PTSD_GSA-MD_v1_3/RP-1194/processamento_inpd2019/proc_inicial_5548inds_dez19/

### IDs estao trocados como enviamos para o Broad. Precisa atuazliar IDs

# Instalar o gatk 

# Criar uma pasta do projeto no seu home

# Entrar na pasta deste novo projeto

# Fazer download e dezipar arquivos do gatk
wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
unzip gatk-4.2.0.0.zip
cd gatk-4.2.0.0/
./gatk 

### Se funcionar o programa vai rodar com as opcoes  

# Pegar dados do cerebro
scp escience@172.23.16.216:/genetica_1/Genotipagem_INPD_Broad2019/PGC_NIMH_Salum_PTSD_GSA-MD_v1_3/RP-1194/processamento_inpd2019/proc_inicial_5548inds_dez19/INPD_GenoBroad19_11.vcf.gz ./IC_JA/
senha: escience

### (scp quem praonde)

# Indexar o vcf. Imagina que esse arquivo é gigantesco, "indexar" em genética significa que o programinha vai criar um index para facilitar encontrar as infos mais rapidamente no arquivo. Tipo um livro
for i in {11..11}; do ./gatk IndexFeatureFile -I /home/julia.arendt/IC_JA/INPD_GenoBroad19_${i}.vcf.gz 
for i in {11..11}; do tabix -p vcf INPD_GenoBroad19_${i}.vcf.gz; done
 
### Cuidado com os caminhos  

# Separar as informações que são importantes para extrair o dado de CNV. 
for i in {11..11}; do ./gatk VariantsToTable -V /home/julia.arendt/IC_JA/INPD_GenoBroad19_${i}.vcf.gz -F ID -F CHROM -F POS -GF BAF -GF LRR -O INPD_CNVs_${i}.vcf; done

# trocando cabeçalhos (isso pode ser feito de uma forma mais bonita)
for i in {11..11}; do sed 's/BAF/B Allele Freq/g' INPD_CNVs_${i}.vcf > INPD_CNVs_${i}.baf ; done
for i in {11..11}; do sed 's/LRR/Log R Ratio/g' INPD_CNVs_${i}.baf > INPD_CNVs_${i}.baf_lrr ; done

## Semelhante ao script da Malu daqui pra baixo. Ver e inserir as hashtags dela aqui

perl ~/PennCNV-1.0.5/compile_pfb.pl 

### Só da primeira vez

for i in {11..11}; do perl ~/IC_JA/PennCNV-1.0.5/kcolumn.pl INPD_CNVs_${i}.baf_lrr split 2 -tab -head 3 -name -out INPD_${i}_ref ; done

# Mudar cabeçalhos dos *_ref* ## Perguntar pro Marquinhos (?)

# Criar pfb unico
for i in {11..11}; do ls INPD_${i}_ref* > IDs_INPD_${i} ; done
#for i in {11..11}; do head -n 250 IDs_INPD_${i} > IDs_ParaPFB ; done 

### Checar um dia se os primeiros 250 não são aparentados. Ver no espelho do Broad.

for i in {11..11}; do sed 's/INPD_${i}_ref.//g' IDs_INPD_${i} > lista_files_INPD_${i} ; done

# mudar cabecalho só do primeiro arquivo, de ID para Name, de POS para Position e de CHROM para Chr
# Para sair do vi salvando: ESC ":x"
# Para sair do vi sem salvar: ESC ":q!"

#vi INPD_11_ref.201904760146_R07C01
#for i in {11..11}; do perl /home/julia.arendt/IC_JA/PennCNV-1.0.5/compile_pfb.pl -list IDs_ParaPFB -output INPD.pfb ; done
# reverter cabecalho para o original
#vi INPD_11_ref.201904760146_R07C01

# Daqui pra baixo você vai preciar mudar os radicais dos arquivos e com os dados da illumina
# Aqui a gente precisava mudar o cabeçalho de todos os arquivos (500 ) mas achamos mais facil mudar a chamada no script, ao inves de procurar por Name, ele procura por ID.
# sed 's/Name/ID/g; s/Chr/CHROM/g; s/Position/POS/g' detect_cnv.pl > detect_cnv_header.pl

for i in {11..11}; do detect_cnv_header.pl -test -hmm /home/julia.arendt/IC_JA/PennCNV-1.0.5/lib/hhall.hmm -pfb /home/julia.arendt/IC_JA/gatk-4.2.0.0/INPD.pfb -list /home/julia.arendt/IC_JA/gatk-4.2.0.0/IDs_INPD11 -log INPD11.log -out INPD11.rawcnv; done

### Mudar sd22 para INPD

### Dentro da pasta gatk

./filter_cnv.pl sd22_affy.rawcnv -qclogfile sd22_affy.log -qclrrsd 0.3 -qcpassout sd22_affy.qcpass -qcsumout sd22_affy.qcsum -qcnumcnv 100 -out sd22_affy.goodcnv



PennCNV - controle de qualidade (Malú)

# 1. QC por amostra
./filter_cnv.pl ../gatk-4.2.0.0/sd22_affy.rawcnv -qclogfile ../gatk-4.2.0.0/sd22_affy.log -qclrrsd 0.3 -qcpassout sd22_affy.qcpass -qcsumout sd22_affy.qcsum -qcnumcnv 100 -out sd22_affy.goodcnv

(wc) 814   5698 120141 ../gatk-4.2.0.0/sd22_affy.rawcnv
(wc) 326  2282 48156 sd22_affy.goodcnv

(wc)  10459   73213 1578334 /home/julia.arendt/IC_JA/gatk-4.2.0.0/INPD11.rawcnv
(wc) 6683   46781 1011934 INPD11.goodcnv

## No R, usando o arquivo sd22_affy.qcsum. Visualização da distribuição de CNVs por pacientes (não fiz!)
eshg <- read.table("~/caminho/arquivoquevcquerchamar.goodcnv", h=T) 
head(eshg)
barplot(eshg$LRR_SD, ylim = c(0, 0.3), xlab = "Pacientes", ylab = "LRR_SD", main = "Sample QC") 
barplot(eshg$NumCNV, ylim = c(0,30), main = "CNVs por amostra", xlab = "Pacientes", ylab = "Número de CNVs" )
0
# 2. Remover regiões ruins

## Centroméricas
### Arquivo centromerescerto.txt achado na internet
# awk '{print $1":"$2"-"$3}' centromere.txt > centromerescerto.txt

./scan_region.pl sd22_affy.goodcnv centromerescerto.txt -minqueryfrac 0.5 > cnvcall.cen
fgrep -v -f cnvcall.cen sd22_affy.goodcnv > sd22_affy_cen.clean

(wc) 322  2254 47557 sd22_affy_cen.clean

(wc) 6638   46466 1004985 INPD11_cen.clean

## Teloméricas:
### Arquivo telomeres.txt feito com 10000bp de diferença - Genome Browser

./scan_region.pl sd22_affy_cen.clean telomeres.txt -minqueryfrac 0.5 > cnvcall.tel
fgrep -v -f cnvcall.tel sd22_affy_cen.clean > sd22_affy_cen_tel.clean

(wc) 320  2240 47266 sd22_affy_cen_tel.clean

(wc) 6638   46466 1004985 INPD11_cen_tel.clean

## Duplicações segmentais, exceto 22q11.2: 
### Arquivo segdup_hg19certo.bed feito em https://genome.ucsc.edu/cgi-bin/hgTables (group: Repeats; tracks: Segmental Dups; table: genomicSuperDups)
tail -n +2 segdup.txt | awk '{print $2,$3, $4}'| sed 's/ /\t/g' > segdup_hg19.bed
awk '{print $1":"$2"-"$3}' segdup_hg19.bed > segdup_hg19certo.bed

./scan_region.pl sd22_affy_cen_tel.clean segdup_hg19certo.bed -minqueryfrac 0.5 > cnvcall.segdup
fgrep -v -f cnvcall.segdup sd22_affy_cen_tel.clean > sd22_affy_cen_tel_segdup.clean

(wc) 258  1806 37974 sd22_affy_cen_tel_segdup.clean

(wc)  5372  37604 811895 INPD11_cen_tel_segdup.clean

## Regiões de imunoglobulinas:
### Não achei arquivo para hg19

./scan_region.pl sd22_affy_cen_tel_segdup.clean immuno_region -minqueryfrac 0.5 > cnvcall.immuno 
fgrep -v -f cnvcall.immuno sd22_affy_cen_tel_segdup.clean > sd22_affy_cen_tel_segdup_immuno.clean

# 3. Merge
./clean_cnv.pl -signalfile ashg.329 combineseg sd22_affy_cen_tel_segdup_immuno.clean > sd22_affy_cen_tel_segdup_immuno_merged.clean

signalfile: a file that contains Chr and Position for SNPs/markers
Observação: fraction 0.2 (default)

### Não fiz

# 4. Mínimo de sondas:

## Deleções (20 sondas):
./filter_cnv.pl -numsnp 20 -length 1k -type del sd22_affy_cen_tel_segdup_immuno_merged.clean -output sd22_affy_del.clean

(wc) 50  350 7307 sd22_affy_del.clean

(wc)  170  1190 25267 INPD11_del.clean

## Duplicações (20 sondas):
./filter_cnv.pl -numsnp 20 -length 1k -type dup sd22_affy_cen_tel_segdup_immuno_merged.clean -output sd22_affy_dup.clean

(wc) 45  315 6684 sd22_affy_dup.clean

(wc)  144  1008 21702 INPD11_dup.clean

# 5. Juntar arquivos com deleções e duplicações:
cat sd22_affy_del.clean sd22_affy_dup.clean > sd22_affy_del_dup.clean

(wc) 95   665 13991 sd22_affy_del_dup.clean

(wc)  314  2198 46969 INPD11_del_dup.clean

# 6. Ordenar pela primeira coluna
man sort 
sort -rk 1 sd22_affy_del_dup.clean > deldup.ordenado





