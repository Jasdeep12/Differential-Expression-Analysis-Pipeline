import click
import subprocess
import re
import sys
from pathlib import Path
import pandas as pd


SAMPLES_TSV = Path('config/samples.tsv')
ACCESSION_PATTERN = re.compile(r'^(SRR|ERR|DRR)\d+$')
FASTQ_SUFFIXES = ('.fastq','.fastq.gz','.fq','.fq.gz')


@click.group()
def cli():
	"""Adding samples for the RNASeq pipeline"""
	pass


@cli.command()
@click.argument('accession', nargs=-1)
@click.option('--name', type=str, default=None, multiple=True, help='Sample names to use in the samples.tsv, enter them in the same order as the accession numbers (defaults to the accession')
@click.option('--force', is_flag=True, help='Overwrite any existing samples with the same name')
def fetchall(accession, name, force):
	"""Downloads multiple ACCESSIONS via prefetch and fasterq-dump and registers it as samples"""
	
	if len(accession) < 2:
		raise click.UsageError(
			"Less than 2 codes detected, for single codes use 'fetch' instead")
	
	if len(name) != len(accession) and len(name) > 0:
		raise click.UsageError(
			"Every sample should have an assigned name.")
		
	
	for code in accession:
		if not ACCESSION_PATTERN.match(code):
			raise click.BadParameter(
			f"'{code}' doesn't seem to be a valid SRA accession"
			f"(expected SRR/ERR/DRR followed by digits)")
	
	raw_dir = Path('data/raw')
	raw_dir.mkdir(parents=True, exist_ok=True)
	failed = []
	
	
	for index, code in enumerate(accession):
		sample_name = name[index] if name else code
		
		try:
			subprocess.run(['prefetch', code, '-O', str(raw_dir)],
				check=True,
				capture_output=True,
				text=True)

		except subprocess.CalledProcessError as e:
			click.echo(f"prefetch failed for {code}", err=True)
			click.echo(f"Exit Code: {e.returncode}", err=True)
			click.echo(f"Error Output / stderr : {e.stderr}")
			click.echo(f"Skipping {code}")
			failed.append(code)
			continue

		except FileNotFoundError:
			click.echo("prefetch not found on PATH. Is the SRA Toolkit activated in your conda env?", err=True)
			sys.exit(1)


		# prefetch will usually drop <accession>.sra inside a subfolder of -O
		sra_path = raw_dir / code / f"{code}.sra"
		if not sra_path.exists():
			sra_path = raw_dir / f"{code}.sra"


		try:
			subprocess.run(['fasterq-dump', str(sra_path), '-O', str(raw_dir)],
				check=True,
				capture_output=True,
				text=True)

		except subprocess.CalledProcessError as e:
			click.echo(f"fasterq-dump failed for {code}", err=True)
			click.echo(f"Exit Code: {e.returncode}", err=True)
			click.echo(f"stderr: {e.stderr}", err=True)
			click.echo(f"Skipping download for {sra_path}")
			failed.append(code)
			continue

		r1 = raw_dir / f"{code}_1.fastq"
		r2 = raw_dir / f"{code}_2.fastq"
		single = raw_dir / f"{code}.fastq"

		if r1.exists() and r2.exists():
			register(sample_name, r1, r2, force=force)
			click.echo(f"Registered '{sample_name}' -> {r1}, {r2}")
		elif single.exists():
			register(sample_name, single, force=force)
			click.echo(f"Registered '{sample_name}' (single-end) -> {single}")
		else:
			click.echo(f"No recognizable fasterq-dump output found for {code}", err=True) 
			failed.append(code)
			continue
			
	if len(failed) > 0:
		message = ', '.join(failed)
	else:
		message = "None"
	
	click.echo(f"{len(accession)-len(failed)}/{len(accession)} processed. Failed: {message}")
		
		


@cli.command()
@click.argument('accession')
@click.option('--name', type=str, default=None, help='Sample name to use in samples.tsv (defaults to the accession)')
@click.option('--force', is_flag=True, help='Overwrite an existing sample with the same name')
def fetch(accession,name,force):
	"""Downloads the ACCESSION via prefetch and fasterq-dump and registers it as a sample"""
	if not ACCESSION_PATTERN.match(accession):
		raise click.BadParameter(
			f"'{accession}' doesn't seem to be a valid SRA accession"
			f"(expected SRR/ERR/DRR followed by digits)")
			

	sample_name = name or accession
	raw_dir = Path('data/raw')
	raw_dir.mkdir(parents=True, exist_ok=True)
	

	try:
		subprocess.run(['prefetch', accession, '-O', str(raw_dir)],
			check=True,
			capture_output=True,
			text=True)

	except subprocess.CalledProcessError as e:
		click.echo(f"prefetch failed for {accession}", err=True)
		click.echo(f"Exit Code: {e.returncode}", err=True)
		click.echo(f"Error Output / stderr : {e.stderr}")
		sys.exit(1)

	except FileNotFoundError:
		click.echo("prefetch not found on PATH. Is the SRA Toolkit activated in your conda env?", err=True)
		sys.exit(1)


	# prefetch will usually drop <accession>.sra inside a subfolder of -O
	sra_path = raw_dir / accession / f"{accession}.sra"
	if not sra_path.exists():
		sra_path = raw_dir / f"{accession}.sra"


	try:
		subprocess.run(['fasterq-dump', str(sra_path), '-O', str(raw_dir)],
			check=True,
			capture_output=True,
			text=True)

	except subprocess.CalledProcessError as e:
		click.echo(f"fasterq-dump failed for {accession}", err=True)
		click.echo(f"Exit Code: {e.returncode}", err=True)
		click.echo(f"stderr: {e.stderr}", err=True)
		sys.exit(1)

	r1 = raw_dir / f"{accession}_1.fastq"
	r2 = raw_dir / f"{accession}_2.fastq"
	single = raw_dir / f"{accession}.fastq"

	if r1.exists() and r2.exists():
		register(sample_name, r1, r2, force=force)
		click.echo(f"Registered '{sample_name}' -> {r1}, {r2}")
	elif single.exists():
		register(sample_name, single, force=force)
		click.echo(f"Registered '{sample_name}' (single-end) -> {single}")
	else:
		click.echo(f"No recognizable fasterq-dump output found for {accession}", err=True)
		sys.exit(1)	 # or append to `failed` + continue, in fetchall's case


@cli.command()
@click.argument('paths', nargs=-1, type=click.Path(exists=True, path_type=Path))
@click.option('--force', is_flag=True, help='Overwrite existing samples with the same name')
def local(paths, force):
	"""Register already downloaded FASTQ files (or directories of them) as samples."""
	if not paths:
		raise click.UsageError("Provide at least one FASTQ file or directory")

	files = []
	for path in paths:
		if path.is_file():
			if path.name.endswith(FASTQ_SUFFIXES):
				files.append(path)
		elif path.is_dir():
			for file in sorted(path.iterdir()):
				if file.is_file() and file.name.endswith(FASTQ_SUFFIXES):
					files.append(file)

	if not files:
		raise click.UsageError("No FASTQ files in the given directories")

	pairs = pairmates(files)


	for sample_name, (r1, r2) in pairs.items():
		if r1 is not None and r2 is not None:
			register(sample_name, r1, r2, force=force)
			click.echo(f"Registered {sample_name}: {r1} / {r2}")
		elif r1 is not None or r2 is not None:
			single = r1 or r2
			register(sample_name, single, force=force)
			click.echo(f"Registered {sample_name} (single-end): {single}")




def pairmates(files):
	"""Group a flat list of FASTQ files into {sample_name: (R1_Path,R2_Path)}
	by detecting _1/_2 or _R1/_R2 in each filename"""
	mate_pattern = re.compile(r'^(.+?)_R?([12])(?:\.fastq|\.fq)(?:\.gz)?$')
	pairs = {}

	for file in files:
		match = mate_pattern.match(file.name)
		if not match:
			click.echo(f"Skipping '{file.name}': couldn't detect R1/R2 in filename", err=True)
			continue
		sample_name, mate = match.group(1), match.group(2)
		pairs.setdefault(sample_name, [None, None])
		pairs[sample_name][0 if mate == '1' else 1] = file

	return {name: tuple(mates) for name, mates in pairs.items()}



def register(sample_name, r1, r2=None, force=None):
	"""Registers files in the data directory with the samples.tsv, allowing them to processed"""

	r1 = Path(r1)
	r2 = Path(r2) if r2 is not None else None

	if not r1.exists():
		raise FileNotFoundError(f"R1 file does not exist: {r1}")


	SAMPLES_TSV.parent.mkdir(parents=True, exist_ok=True)

	if SAMPLES_TSV.exists():
		df = pd.read_csv(SAMPLES_TSV, sep='\t')
	else:
		df = pd.DataFrame(columns=['sample','R1','R2'])

	r2_value = str(r2) if r2 is not None else ''
	existing = df['sample'] == sample_name

	if existing.any():
		if not force:
			raise ValueError(f"Sample: {sample_name} already exists in samples.tsv. Use --force to Overwrite.")
		df.loc[existing, ['R1','R2']] = [str(r1), r2_value]
	else:
		new_row = pd.DataFrame([{'sample': sample_name, 'R1': str(r1), 'R2': r2_value}])
		df = pd.concat([df, new_row], ignore_index=True)

	df.to_csv(SAMPLES_TSV, sep='\t', index=False)


@cli.command()
@click.option('--control',type=int, default=0,help='Number of control samples (Assuming those are the first N entries)')
@click.option('--treatment', multiple=True, type=int, help='Number of treatment samples, can have multiples samples (Assuming they are in correct order), Ex. Entries 4-6 are treatment 1, and Entries 7-9 are treament 2, So "--treatment 3 --treatment 3"')
@click.option('--names', multiple=True, type=str, help='Names of treatment(s), Ex. "--names treatment1 --names treatment2"')
def condition(control,treatment,names):
	"""Add conditions to the samples.tsv for each respective entry."""

	df = pd.read_csv(SAMPLES_TSV, sep='\t')
	if df.empty:
		raise click.UsageError(f"{SAMPLES_TSV} is empty, please add samples.")

	if len(df) != (control + sum(treatment)):
		raise click.UsageError(f"Number of control samples ({control}) and number of treatment samples ({sum(treatment)}) does not match number of samples found in {SAMPLES_TSV}")

	if len(names) != len(treatment):
		raise click.BadParameter(f"Number of treatments ({len(treatment)}) does not match number of names given ({len(names)})")

	columns = []
	for count in range(control):
		columns.append('control')
	for count, name in zip(treatment, names):
		for j in range(count):
			columns.append(f'{name}')

	df['condition'] = columns
	click.echo(f'Condition assignment for {SAMPLES_TSV} is as follows:')
	click.echo(df[['sample','condition']].to_string(index=False))

	df.to_csv(SAMPLES_TSV, sep='\t',index=False)



@cli.command()
@click.option('--samples',is_flag=True,default=False,help=f'Clears {SAMPLES_TSV}')
@click.option('--raw-data',is_flag=True,default=False,help=f'Clears data/raw')
@click.option('--results',is_flag=True,default=False,help=f'Clears results/')
@click.option('--logs',is_flag=True,default=False,help=f'Clears Logs/')
def clear(samples,raw_data,results,logs):
	"""Easy way to clean the pipeline for reuse, or debugging"""
	pass
	


if __name__ == '__main__':
	cli()












