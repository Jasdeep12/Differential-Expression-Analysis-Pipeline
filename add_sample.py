import click
import subprocess
from pathlib import Path

@click.group()

def cli():
	"""Adding samples for the RNASeq pipeline"""
	pass



@cli.command()
@click.option('--fetch', type=str, help='Accession Number')
def fetchAccession(fetch):
	"""Fetches the accession number and adds it data directory"""

	if fetch != str:
		raise ValueError(f"{fetch} is not a valid accession number")

	try:
		accession = subprocess.run(['prefetch', fetch, '-O', 'data/raw/'],
			check=True,
			capture_output=True)

	except subprocess.CalledProcessError() as e:
		print("Fetching Failed")
		print(f"Exit Code: {e.returncode}")
		print(f"Error Output: {e.stderr}")

	try:
		subprocess.run(['fasterq-dump', f"data/raw/{fetch}.sra"],
			check=True,
			capture_output=True)

	except subprocess.CalledProcessError() as e:
		print("fasterq-dump failed")
		print(e.stderr)

	#Todo Build file register for Sample.tsv

@cli.command()
@click.option('--local', type=click.Path(exists=True,file_okay=True, dir_okay=True, path_type=Path), multiple=True, help='Add FASTQ files')
def addFile(local):
	if local == None:
		raise FileNotFoundError

	for f in local:


