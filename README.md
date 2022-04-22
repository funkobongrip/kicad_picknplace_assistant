# KiCAD PickNPlace Assistant

Manual pick-and-place document generator for KiCAD PCB files.

A simple python script that generates PDF files that indicate the position of components on the PCB for every item on the BOM list.
Each BOM line generates one corresponding page.

For a more advanced, interactive solution see: [InteractiveHtmlBom](https://github.com/openscopeproject/InteractiveHtmlBom)

## Getting Started and Usage

Simply call the script with Python 3, and provide the filename of your *.kicad_pcb file.
```
python kicad_picknplace_assistant.py your_file.kicad_pcb
```

Added a few features:

  -s, --split           split into one file per component

  -p, --precise         be precise with footprints (else C0402N == C0402L)

  -m, --dontconnectorgroup ignore connector values when grouping parts

  -f, --folders         sort pdfs (when splitting) to folders tht smt smt/passives

  -n BOARDS, --panels=BOARDS set number of boards/panels to populate

  -b, --bom-only        print only bom, then exit

  -c, --csv-bom         create csv bom file


### Prerequisites

This script runs on Python 3.

Requires matplotlib, numpy etc.

Requires the KiCAD's pcbnew to be installed (it comes with the required Python libraries).

## Acknowledgments

Original script by [@pwuertz](https://github.com/pwuertz)
