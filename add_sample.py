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
@click.argument('accession')
@click.option('--name', type=str, default=None, help='Sample name to use in samples.tsv (defaults to the accession)')
@click.option('--force', is_flag=True, help='Overwrite an existing sample with the same name')
def fetch(accession,name,force):
	"""Downloads the ACCESSION via prefetch and fasterq-dump and registers it as a sample"""
	if not ACCESSION_PATTERN.match(accession):
		raise click.BadParameter(
			f"'{accession}' doesn't seem to be a valid SRA accession")
			f"(expected SRR/ERR/DRR followed by digits)"
			

	sample_name = name or accession
	raw_dir = Path('data/raw')
	raw_dir.mkdir(parents=True, exist_ok=True)
	

	try:
		accession = subprocess.run(['prefetch', fetch, '-O', str(raw_dir)],
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

	if not r1.exists() or not r2.exists():
		click.echo(f"Expected paired output not found: {r1} / {r2}", err=True)
		click.echo("Check whether this is single-end data, or fasterq-dump used a different naming pattern.", err=True)
		sys.exit(1)

	register(sample_name, r1, r2, force=force)
	click.echo(f"Registered '{sample_name}' -> {r1}, {r2}")


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
			if path.endswith(FASTQ_SUFFIXES):
				files.append(path)
		elif path.is_dir():
			for file in sorted(path.iterdir()):
				if file.is_file and file.endswith(FASTQ_SUFFIXES):
					files.append(file)

	if not files:
		raise click.UsageError("No FASTQ files in the given directories")

	pairs = pairmates(files)


	for sample_name, [r1,r2] in pairs.items():
		if r1 is None or r2 is None:
			missing_mate = 'R2' if r2 is None or 'R1' if r1 is None
			click.echo(f"{sample_name} is missing {missing_mate}: skipping {sample_name}", err=True)
			continue
		register(sample_name, r1, r2)
		click.echo(f"Registered {sample_name} : {r1} / {r2} ")





def pairmates(files):
	"""Group a flat list of FASTQ files into {sample_name: (R1_Path,R2_Path}
	by detecting _1/_2 or _R1/_R2 in each filename"""
	mate_pattern = re.compile(r'^(.+?)_R?([12])(?:\.fastq|\.fq)(?:\.gz)?$')
	pairs = {}

	for file in files:
		match = mate_pattern.match(file.name)
		if not match:
			click.echo(f"Skipping '{file.name}: couldn't detect R1/R2 in filename", err=True)
			continue
		sample_name, mate = match.group(1), match.group(2)
		pairs.setdefault(sample_name, [None,None])
		pairs[sample_name][0 if mate == '1' else 1] = file

	return {name: tuple(mates) for name, mates in pair.items()}



def register(r1, r2, sample_name):
	"""Registers files in the data directory with the samples.tsv, allowing them to processed"""
	samples = Path('config/samples.tsv')

	r1, r2 = Path(r1), Path(r2)

	if not r1.exists():
		raise FileNotFoundError(f"R1 file does not exist: {r1}")

	if not r2.exists():
		raise FileNotFoundError(f"R2 file does not exist: {r2}")

	SAMPLES_TSV.parent.mkdir(parent=True, exist_ok=True)

	if SAMPLES_TSV.exists():
		df = pd.read_csv(SAMPLES_TSV, sep='\t')
	else:
		df = pd.DataFrame(columns=['sample','r1','r2'])

	existing = df['sample'] == sample_name

	if existing.any():
		if not force:
			raise ValueError(f"Sample: {sample_name} already exists in samples.tsv. Use --force to Overwrite.")

		df.loc[exising, ['R1','R2']] = [str(r1), str(2)]
	else:
		new_row = pd.DataFrame([{'sample': sample_name, 'R1': str(r1), 'R2': str(r2)}])
		df = pd.concat([df,new_row], ignore_index=True)

	df.to_csv(SAMPLES_TSV, sep='\t', index=False)


if __name__ == '__main__':
	cli()















