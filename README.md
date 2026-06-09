# ANNIE Helper Scripts

A repository for small standalone scripts that are detached from the ToolAnalysis. 
For debugging, validation, ML training, and similar scripts. 

Please note that ToolAnalysis remains the primary software framework for ANNIE. This repository is intended to host standalone helper scripts and utilities that complement the main codebase.

## Repository Structure

```text
.
├── cpp/
├── docs/
│   ├── script_catalog.md
│   └── script_header_template.md
├── python/
├── .gitignore
└── README.md
```

## Directory Description

| Directory | Description                                             |
| --------- | ------------------------------------------------------- |
| `cpp/`    | C++ analysis and utility scripts                        |
| `python/` | Python analysis and utility scripts                     |
| `docs/`   | Repository documentation, templates, and script catalog |

## Documentation Standards

Each script directory should contain:

* The script source code.
* A dedicated `README.md`.
* A standardized script header.

Templates are available in:

* `docs/script_header_template.md`

## Script Catalog

A catalog of all documented scripts can be found in:

* `docs/script_catalog.md`

## Repository Guidelines

When adding new scripts to this repository:

* Place scripts in the appropriate sub-directory (`cpp/` or `python/`).
* The sub-directory name should reflect the primary purpose of the script(s) it contains.
* All scripts must include the standard metadata header defined in `docs/script_header_template.md`.
* Each script directory must contain a `README.md` describing the script, its inputs, outputs, dependencies (if applicable), and usage.
* Add all new scripts to `docs/script_catalog.md`.
* Update documentation whenever script functionality changes.
* For c++ scripts, please include a Makefile to compile your application, unless it is to run as ROOT macro (though compiled scripts are preferred for compile-time error checking).
* For python scripts, please specify your environment and all python modules, including versions which your script relies on.
