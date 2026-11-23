# intermediate-tail-change-date

R code and synthetic simulation materials for *Estimated marginal change dates and the size of an intermediate-tail dependence test* by Çağlar SÖZEN and Mervenur SÖZEN.

The 24,000 evaluation rows and 456,000 associated calibration statistics come from the fixed design in the article. Each evaluation uses the same cutoff for the three date rules. For a nominal level of 0.05, that cutoff is the maximum of its 19 fresh true-date calibration statistics. All observations were generated from the stated Clayton model.

The exact seed streams from the completed run are in [data/seeds.rds](data/seeds.rds). The file [data/calibration_values.csv.gz](data/calibration_values.csv.gz) contains all calibration statistics with identifiers for design cell, replication, and calibration draw. The raw uniforms can be regenerated from the seed states. [data/trials.csv](data/trials.csv) contains the evaluation decisions and diagnostics. The original design and versions are recorded in [FULL_PROTOCOL.json](FULL_PROTOCOL.json), [data/executed_design.R](data/executed_design.R) and [data/session_info.txt](data/session_info.txt). [provenance.json](provenance.json) records the file digests and the original run identifier.

Use the repository root as the working directory. The archived run used base R 4.5.0. To regenerate one original evaluation and its 19 calibration statistics from their assigned streams, run

    Rscript --vanilla reproduce.R chunk 1 1 1 replay_one

The three numbers are the cell, first replication and count. Each call handles at most ten replications and writes to a new folder. To regenerate the summary tables directly from the recorded evaluation rows, run

    Rscript --vanilla reproduce.R tables rebuilt_tables

The support column in the numerical files uses the known zero endpoint. The prefix extension in the article is a theoretical result.

See [LICENSE](LICENSE) for the code license.
