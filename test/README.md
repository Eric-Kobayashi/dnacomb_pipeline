# Test data
Simple test data to run all stages of the pipeline on. Useful for testing an environment
is set up correctly or when developing new builds.

Run the test with Singularity:

`make test`

To choose a different container runtime, run Nextflow directly with
`-config test/test.config`.
