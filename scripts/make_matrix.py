import pandas as pd
import sys



def make_matrix(inputFile,sampleNames,outputFile):
	merged = None
	for filepath, name in zip(inputFiles, sampleNames):

		df = pd.read_csv(filepath, sep='\t', comment='#')

		countCol = df.columns[-1]

		df = df[['Geneid',countCol]].rename(columns = {countCol: name})

		if merged == None:
			merged = df
		else:
			merged = merged.merge(df, on='Geneid',how='outer')

	merged.to_csv(outputFile,sep='\t',index=False)

if __name__ == '__main__':
	make_matrix(
		inputFile=snakemake.input.counts,
		sampleNames=snakemake.params.samples,
		outputFile=snakemake.output[0]


		)
