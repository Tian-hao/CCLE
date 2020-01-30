#!/usr/bin/env python
from scipy.stats import ranksums
from scipy.stats import fisher_exact
from scipy.stats.mstats import gmean


def main():
  sumfile = open('summary_CCLE_fisher_test.txt','w')
  testlist = ['both','STK11','KRAS']
  testlist = ['STK11','KRAS']
  for test_group in testlist:
    prepare_hetmap('5gene','noSCC',sumfile,test_group)
    prepare_hetmap('5gene','Adeno',sumfile,test_group)
    prepare_hetmap('7gene','noSCC',sumfile,test_group)
    prepare_hetmap('7gene','Adeno',sumfile,test_group)
  sumfile.close()
  

def prepare_hetmap(gene_annot,cell_annot,sumfile,test_group):
  reffile = open('Cell_lines_annotations_20181226.txt')
  header = reffile.readline().rstrip().rsplit('\t')
  hindex = header.index('PATHOLOGIST_ANNOTATION')
  celllist = []
  for line in reffile:
    line = line.rstrip().rsplit('\t')
    if cell_annot == 'Adeno':
      if 'LUNG' in line[0] and 'Adeno' in line[hindex]:
        celllist.append(line[0])
    if cell_annot == 'noSCC':
      if 'LUNG' in line[0] and 'Mesothelioma' not in line[hindex] and 'Lung:SCLC' not in line[hindex]:
        celllist.append(line[0])
  reffile.close()
  #if 'A549_LUNG' not in celllist:
  #  celllist.append('A549_LUNG')

  mutfile = open('CCLE_DepMap_18q3_maf_20180718.txt')
  mutdict = {}
  mutlist = ['KRAS','STK11','both','either']
  header  = mutfile.readline().rstrip().rsplit('\t')
  colcount = header.index('Tumor_Sample_Barcode')
  for line in mutfile:
    line = line.rstrip().rsplit('\t')
    cell = line[colcount]
    if 'LUNG' not in cell: continue
    if cell not in mutdict: mutdict[cell] = {}
    if line[0] in mutlist: mutdict[cell][line[0]] = 1
  mutfile.close()
  mutdict['NCIH1568_LUNG']['STK11'] = 1
  mutdict['NCIH2126_LUNG']['STK11'] = 1
  for cell in celllist:
    if cell not in mutdict: mutdict[cell] = {}
    mutdict[cell]['both'] = 1; mutdict[cell]['either'] = 0
    for mut in mutlist[:2]:
      if mut not in mutdict[cell]:
        mutdict[cell][mut] = 0; mutdict[cell]['both'] = 0
      else:
        mutdict[cell]['either'] = 1  

  rnafile = open('CCLE_RNAseq_genes_rpkm_20180929.gct')
  if gene_annot == '7gene':
    genelist = ['CXCL1','CXCL2','CXCL3','CXCL5','CXCL6','PPBP','IL8']
  if gene_annot == '5gene':
    genelist = ['CXCL1','CXCL2','CXCL3','CXCL5','IL8']
  #genelist = ['CSF1','CSF2','CSF3']
  rnadict = {}
  header = rnafile.readline();header = rnafile.readline();header = rnafile.readline().rstrip().rsplit('\t')
  liblist = header[2:]
  celllist = [x for x in celllist if x in liblist]
  for cell in celllist:
    rnadict[cell] = {}
  for line in rnafile:
    line = line.rstrip().rsplit('\t')
    if line[1] not in genelist: continue
    for i,count in enumerate(line[2:]):
      if liblist[i] in celllist:
        rnadict[liblist[i]][line[1]] = float(count)
  rnafile.close()
  #print(rnadict[next(iter(rnadict))])

  outfile = open('CCLE_combine_'+gene_annot+'_'+cell_annot+'.txt','w')
  outfile.write('cells\t'+'\t'.join(celllist)+'\n')
  for gene in genelist:
    outfile.write(gene)
    for cell in celllist:
      outfile.write('\t'+str(rnadict[cell][gene]))
    outfile.write('\n')
  for mut in mutlist[:2]:
    outfile.write(mut)
    for cell in celllist:
      if cell not in mutdict or mutdict[cell][mut] == 0:
        color = '#808080'
      else:
        color = '#FF0000'
      if cell not in mutdict:
        mutdict[cell] = {}
        for mut in mutlist: mutdict[cell][mut] = 0
      outfile.write('\t'+color)
    outfile.write('\n')
  outfile.close()
  
  #for ranksum test on geneexpression
  ranksumfile = open('CCLE_ranksum_expression'+gene_annot+'_'+cell_annot+'_'+test_group+'.txt','w')
  norm_rnadict = {}; mean_rnadict = {}
  for gene in genelist:
    mean_rnadict[gene] = []
  for cell in rnadict:
    for gene in rnadict[cell]:
      mean_rnadict[gene].append(rnadict[cell][gene]+0.0001)
  for gene in genelist:
    mean_rnadict[gene] = gmean(mean_rnadict[gene])
    #mean_rnadict[gene] = sum(mean_rnadict[gene])/len(celllist)
  for cell in rnadict:
    norm_rnadict[cell] = {}
    for gene in rnadict[cell]:
      norm_rnadict[cell][gene] = rnadict[cell][gene]/mean_rnadict[gene]
  mutrnalist = {}; otherrnalist = {}
  for cell in celllist:
    if mutdict[cell]['both'] == 1: 
      if cell not in mutrnalist: mutrnalist[cell] = 0
      for gene in genelist:
        mutrnalist[cell] += norm_rnadict[cell][gene]
    if mutdict[cell][test_group] == 1 and mutdict[cell]['both'] == 0:
      if cell not in otherrnalist: otherrnalist[cell] = 0
      for gene in genelist:
        otherrnalist[cell] += norm_rnadict[cell][gene]
  ranksumfile.write('cell\tmutation\tsummed_expression\t'+'\t'.join([x+'_normalized' for x in genelist])+'\t'+'\t'.join(genelist)+'\n')
  for cell in mutrnalist:
    ranksumfile.write(cell+'\tdouble\t'+str(mutrnalist[cell]))
    for gene in genelist:
      ranksumfile.write('\t'+str(norm_rnadict[cell][gene]))
    for gene in genelist:
      ranksumfile.write('\t'+str(rnadict[cell][gene]))
    ranksumfile.write('\n')
  for cell in otherrnalist:
    ranksumfile.write(cell+'\t'+test_group+'\t'+str(otherrnalist[cell]))
    for gene in genelist:
      ranksumfile.write('\t'+str(norm_rnadict[cell][gene]))
    for gene in genelist:
      ranksumfile.write('\t'+str(rnadict[cell][gene]))
    ranksumfile.write('\n')
  ranksumfile.close() 
  print(gene_annot,cell_annot,test_group,len(mutrnalist),len(otherrnalist),ranksums(list(mutrnalist.values()),list(otherrnalist.values())))


  #for fisher exact test on clustered data
  dmethodlist = ["euclidean", "maximum", "manhattan", "canberra" ,"minkowski"]
  hmethodlist = ["single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid" ,"median"]
  for dmethod in dmethodlist:
    for hmethod in hmethodlist:
      cxlow = open('CCLE_CXCL_low_'+gene_annot+'_'+cell_annot+'_'+dmethod+'_'+hmethod+'.txt')
      cxhig = open('CCLE_CXCL_high_'+gene_annot+'_'+cell_annot+'_'+dmethod+'_'+hmethod+'.txt')
      oddmat = [[0,0],[0,0]]
      header = cxlow.readline()
      for line in cxlow:
        line = line.rstrip().rsplit('\t')
        cell = line[1]
        if mutdict[cell][test_group] == 1:
          oddmat[0][0] += 1
        else:
          oddmat[0][1] += 1
      header = cxhig.readline()
      for line in cxhig:
        line = line.rstrip().rsplit('\t')
        cell = line[1]
        if mutdict[cell][test_group] == 1:
          oddmat[1][0] += 1
        else:
          oddmat[1][1] += 1
      cxlow.close()
      cxhig.close()
      #print(gene_annot,cell_annot,dmethod,hmethod,oddmat,fisher_exact(oddmat))
      sumfile.write(gene_annot+'\t'+cell_annot+'\t'+dmethod+'\t'+hmethod+'\t'+test_group+'\t'+str(oddmat[0][0])+','+str(oddmat[0][1])+','+str(oddmat[1][0])+','+str(oddmat[1][1])+'\t'+str(fisher_exact(oddmat)[1])+'\n')

if __name__ == '__main__':
  main()
