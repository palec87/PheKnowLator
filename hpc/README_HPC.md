# PheKnowLator Data Preparation - HPC Version

This directory contains the HPC-compatible version of the PheKnowLator data preparation notebook, converted to a Python script with comprehensive logging capabilities.

## Files

- **`data_preparation_hpc.py`** - Main Python script for data preparation
- **`submit_data_prep.sh`** - SLURM batch submission script
- **`README_HPC.md`** - This documentation file

## Features

### ✅ HPC-Ready
- Designed to run on HPC clusters with SLURM scheduling
- Configurable resource requirements
- Background execution support

### ✅ Comprehensive Logging
- Timestamped log files for each run
- Both file and console output
- Progress tracking with tqdm progress bars
- Error handling with full stack traces

### ✅ Flexible Execution
- Run all steps or individual steps
- Command-line argument parsing
- Skip downloads if files already exist

### ✅ Modular Design
- Object-oriented structure
- Easy to extend and maintain
- Reusable components

## Usage

### Local Testing

```bash
# Run all steps
python data_preparation_hpc.py

# Run specific step
python data_preparation_hpc.py --step genomic_ids

# Custom configuration
python data_preparation_hpc.py \
    --log-dir ./custom_logs \
    --data-dir /path/to/data \
    --skip-downloads
```

### HPC Submission (SLURM)

```bash
# Run all steps
sbatch submit_data_prep.sh

# Run specific step
sbatch submit_data_prep.sh genomic_ids

# Check job status
squeue -u $USER

# View log files
tail -f logs/data_prep_<job_id>.out
tail -f logs/data_prep_<job_id>.err
```

## Available Steps

- **`all`** - Run complete data preparation pipeline (default)
- **`genomic_ids`** - Process genomic identifier mappings
- **`mesh_chebi`** - Create MeSH-ChEBI mappings
- **`disease_ids`** - Process disease and phenotype identifiers (not yet fully implemented)
- **`hpa_uberon`** - Process HPA tissue/cell mappings (not yet fully implemented)
- **`reactome_pw`** - Process Reactome pathway mappings (not yet fully implemented)
- **`genomic_so`** - Process genomic sequence ontology mappings (not yet fully implemented)
- **`ontologies`** - Process ontology data (not yet fully implemented)
- **`metadata`** - Generate metadata (not yet fully implemented)

## Command-Line Options

### `--log-dir` (default: `./logs`)
Directory where log files will be written. A timestamped log file will be created for each run.

### `--data-dir` (default: `../resources`)
Base directory for data files. The script expects the following subdirectories:
- `processed_data/unprocessed_data/` - Raw downloads
- `processed_data/` - Processed outputs
- `relations_data/` - Relations data
- `node_data/` - Node metadata
- `construction_approach/` - Construction approach files
- `ontologies/` - Ontology files

### `--step` (default: `all`)
Which processing step to run. See "Available Steps" above.

### `--skip-downloads`
Skip downloading files if they already exist locally. Useful for re-running steps after failures.

## Log Files

Log files are created in the specified log directory with the format:
```
logs/data_preparation_YYYYMMDD_HHMMSS.log
```

Each log entry includes:
- Timestamp
- Module name
- Log level (INFO, WARNING, ERROR)
- Message

Example log output:
```
2026-01-09 14:23:45 - __main__ - INFO - Logging initialized. Log file: logs/data_preparation_20260109_142345.log
2026-01-09 14:23:45 - __main__ - INFO - Script started at 2026-01-09 14:23:45.123456
2026-01-09 14:23:45 - __main__ - INFO - ================================================================================
2026-01-09 14:23:45 - __main__ - INFO - STEP 1: Processing Genomic Identifiers
2026-01-09 14:23:45 - __main__ - INFO - ================================================================================
2026-01-09 14:23:46 - __main__ - INFO - Loading genomic typing dictionary...
...
```

## SLURM Configuration

The `submit_data_prep.sh` script is pre-configured with reasonable defaults:

- **Time limit**: 48 hours
- **Memory**: 64GB
- **CPUs**: 8 cores
- **Partition**: standard (adjust for your cluster)

### Adjusting Resources

Edit the `#SBATCH` directives in `submit_data_prep.sh`:

```bash
#SBATCH --time=72:00:00      # Increase time limit
#SBATCH --mem=128G           # Increase memory
#SBATCH --cpus-per-task=16   # More CPUs
#SBATCH --partition=highmem  # Different partition
```

### Loading Modules

Uncomment and adjust the module loading section for your cluster:

```bash
module load python/3.8
module load conda
source activate ontology
```

## Dependencies

The script requires the same dependencies as the original notebook:

- pandas
- numpy
- networkx
- rdflib
- tqdm
- openpyxl
- reactome2py
- pkt_kg (PheKnowLator package)

Install with:
```bash
pip install -r notebooks/requirements.txt
```

## Progress Tracking

The script provides real-time progress information:

1. **Console output**: Live updates streamed to stdout
2. **Log files**: Persistent record of all operations
3. **Progress bars**: Visual feedback for long-running operations
4. **SLURM logs**: Separate output (.out) and error (.err) files

## Error Handling

- All exceptions are logged with full stack traces
- Script exits with appropriate error codes
- Partial progress is preserved (processed files remain)
- Can resume from where it failed by re-running

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View live output
tail -f logs/data_prep_<job_id>.out

# View errors
tail -f logs/data_prep_<job_id>.err

# Check resource usage
seff <job_id>

# Cancel a job
scancel <job_id>
```

## Performance Considerations

### Memory Usage
- Peak memory usage can exceed 50GB during master dictionary creation
- Consider requesting 64GB+ for safety
- Monitor with `seff <job_id>` after completion

### Time Estimates
- Full pipeline: 8-12 hours
- Genomic identifiers: 4-6 hours (including master dict creation)
- Other individual steps: 30 minutes - 2 hours each

### Optimization Tips
1. Use `--skip-downloads` for reruns to save time
2. Run individual steps in parallel on different nodes
3. Increase CPU count for better tqdm performance
4. Use fast local storage if available

## Troubleshooting

### Job Fails with Out of Memory
- Increase `--mem` in SLURM script
- Check if node has swap space available
- Consider using a high-memory partition

### Download Timeouts
- FTP servers may be slow or temporarily unavailable
- Script will retry automatically
- Check network connectivity from compute nodes

### Import Errors
- Ensure conda environment is activated
- Verify all dependencies are installed
- Check Python version compatibility (3.8+)

### Permission Errors
- Ensure output directories are writable
- Check file permissions for OWLTools
- Verify data directory paths are correct

## Next Steps

After data preparation completes:

1. Verify all output files in `resources/processed_data/`
2. Check log files for any warnings
3. Proceed with knowledge graph construction
4. Use mapping files for edge list generation

## Support

For issues specific to:
- **PheKnowLator**: See [GitHub Issues](https://github.com/callahantiff/PheKnowLator/issues)
- **HPC script**: Create an issue describing the problem with log excerpts
- **SLURM**: Consult your cluster's documentation or support team

## Version History

- **v1.0** (2026-01-09): Initial HPC conversion
  - Implemented genomic identifiers and MeSH-ChEBI steps
  - Added comprehensive logging
  - Created SLURM submission script

## License

Same as PheKnowLator project.
