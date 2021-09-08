

# Os arquivos que você vai usar estão nessa pasta abaixo no Cerebro. São doze arquivos INPD_GenoBroad19_* . Cada um pesa 6.1Gb zipado. Já transferi o primeiro para a sua pasta

/genetica_1/Genotipagem_INPD_Broad2019/PGC_NIMH_Salum_PTSD_GSA-MD_v1_3/RP-1194/processamento_inpd2019/proc_inicial_5548inds_dez19/

# IDs estao trocados como enviamos para o Broad. Precisa atuazliar IDs

# instalar o gatk 

## Criar uma pasta do projeto no seu home

## entrar na pasta deste novo projeto

# Fazer download e dezipar arquivos do gatk
wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip

unzip gatk-4.2.0.0.zip
cd gatk-4.2.0.0/
./gatk # Se funcionar o programa vai rodar com as opcoes  

# indexar o vcf. Imagina que esse arquivo é gigantesco, "indexar" em genética significa que o programinha vai criar um index para facilitar encontrar as infos mais rapidamente no arquivo. Tipo um livro

./gatk IndexFeatureFile -I /home/julia.arendt/IC_JA/INPD_GenoBroad19_12.vcf.gz # cuidado com os caminhos  

# Separa as informações que são importantes para extrair o dado de CNV. Você pode ( e deve) estudar o que signifca cada uma dessas colunas que eu selecionei. Principalmente o BAF e LRR.

./gatk VariantsToTable -V /home/julia.arendt/IC_JA/INPD_GenoBroad19_12.vcf.gz -F ID -F CHROM -F POS -GF BAF -GF LRR -O INPD_CNVs_12.vcf  


# trocando cabeçalhos (isso pode ser feito de uma forma mais bonita)

sed 's/BAF/B Allele Freq/g' INPD_CNVs_12.vcf > INPD_CNVs_12.baf

sed 's/LRR/Log R Ratio/g' INPD_CNVs_12.baf > INPD_CNVs_12.baf_lrr


# Semelhante ao script da Malu daqui pra baixo. Ver e inserir as hashtags dela aqui

perl ~/PennCNV-1.0.5/compile_pfb.pl

perl ~/PennCNV-1.0.5/kcolumn.pl INPD_CNVs_12.baf_lrr split 2 -tab -head 3 -name -out INPD_ref # parar aqui para apróxima semana


ls INPD_ref* > IDs_INPD12

sed 's/INPD_ref.//g' IDs_INPD12 > lista_files_INPD12

perl /home/julia.arendt/IC_JA/PennCNV-1.0.5/compile_pfb.pl -list IDs_INPD12 -output INPD12.pfb


# daqui pra baixo você vai preciar mudar os radicais dos arquivos e com os dados da illumina

detect_cnv.pl -test -hmm /home/julia.arendt/IC_JA/PennCNV-1.0.5/affy/libgw6/affygw6.hmm -pfb /home/julia.arendt/IC_JA/gatk-4.2.0.0/INPD12.pfb -list /home/julia.arendt/IC_JA/gatk-4.2.0.0/IDs_INPD12 -log sd22_affy.log -out sd22_affy.rawcnv

./filter_cnv.pl sd22_affy.rawcnv -qclogfile sd22_affy.log -qclrrsd 0.3 -qcpassout sd22_affy.qcpass -qcsumout sd22_affy.qcsum -qcnumcnv 100 -out sd22_affy.goodcnv

PennCNV - controle de qualidade (Malú)

# 1. QC por amostra
./filter_cnv.pl ../gatk-4.2.0.0/sd22_affy.rawcnv -qclogfile ../gatk-4.2.0.0/sd22_affy.log -qclrrsd 0.3 -qcpassout sd22_affy.qcpass -qcsumout sd22_affy.qcsum -qcnumcnv 100 -out sd22_affy.goodcnv

(wc) 814   5698 120141 ../gatk-4.2.0.0/sd22_affy.rawcnv
(wc) 326  2282 48156 sd22_affy.goodcnv

No R, usando o arquivo sd22_affy.qcsum. Visualização da distribuição de CNVs por pacientes
eshg <- sd22_affy barplot(eshg$LRR_SD, ylim = c(0, 0.3), xlab = "Pacientes", ylab = "LRR_SD", main = "Sample QC") barplot(eshg$NumCNV, ylim = c(0,30), main = "CNVs por amostra", xlab = "Pacientes", ylab = "Número de CNVs" )

# 2. Remover regiões ruins

## Centroméricas
./scan_region.pl sd22_affy.goodcnv centromerescerto.txt -minqueryfrac 0.5 > cnvcall.cen
fgrep -v -f cnvcall.cen sd22_affy.goodcnv > sd22_affy_cen.clean

(wc) 322  2254 47557 sd22_affy_cen.clean

## Teloméricas:
./scan_region.pl sd22_affy_cen.clean telomeres.txt -minqueryfrac 0.5 > cnvcall.tel
fgrep -v -f cnvcall.tel sd22_affy_cen.clean > sd22_affy_cen_tel.clean

(wc) 320  2240 47266 sd22_affy_cen_tel.clean

## Duplicações segmentais, exceto 22q11.2 
./scan_region.pl sd22_affy_cen_tel.clean segdup_hg19certo.bed -minqueryfrac 0.5 > cnvcall.segdup
fgrep -v -f cnvcall.segdup sd22_affy_cen_tel.clean > sd22_affy_cen_tel_segdup.clean

(wc) 258  1806 37974 sd22_affy_cen_tel_segdup.clean

## Regiões de imunoglobulinas:
./scan_region.pl sd22_affy_cen_tel_segdup.clean immuno_region -minqueryfrac 0.5 > cnvcall.immuno 
fgrep -v -f cnvcall.immuno sd22_affy_cen_tel_segdup.clean > sd22_affy_cen_tel_segdup_immuno.clean


# 3. Merge
./clean_cnv.pl -signalfile ashg.329 combineseg sd22_affy_cen_tel_segdup_immuno.clean > sd22_affy_cen_tel_segdup_immuno_merged.clean

signalfile: a file that contains Chr and Position for
SNPs/markers
Observação: fraction 0.2 (default)

# 4. Mínimo de sondas:
Deleções: 20 sondas
./filter_cnv.pl -numsnp 10 -type del sd22_affy_cen_tel_segdup_immuno_merged.clean -output sd22_affy_del.clean

Duplicações: 20 sondas
./filter_cnv.pl -numsnp 20 -type dup sd22_affy_cen_tel_segdup_immuno_merged.clean -output sd22_affy_dup.clean

# 5. Juntar arquivos com deleções e duplicações:
cat sd22_affy_del.clean sd22_affy_dup.clean > sd22_affy_del_dup.clean

# 6.Ordenar pela primeira coluna
man sort 
sort -rk 1 nomedoarquivo > nomedoarquivonovo_ordenado





