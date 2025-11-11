# Humatch

**Fast, gene-specific joint humanisation of antibody heavy and light chains.**

```
                        @
     )  __QQ    -->    /||\
    (__(_)_">           /\         
```

<!--- INSTALL --->
## Install

Follow the steps below to download the code and the necessary packages. Requires **Python 3.9**.

```
# clone the repo
git clone https://github.com/oxpig/Humatch.git
cd Humatch/

# create your virtual env e.g.
python3 -m venv .humatch_venv
source .humatch_venv/bin/activate

# install
pip install .
```

**ANARCI** is required for aligning and padding sequences. We recommend installing ANARCI from <a href="https://github.com/oxpig/ANARCI/tree/master">github.com/oxpig/ANARCI</a>.

If you are having issues installing, try upgrading pip: ```pip install --upgrade pip```

Specific python versions can be used to initiate the environment using e.g. ```/usr/bin/python3.9 -m venv .humatch_venv```

CNN weights and germline likeness lookup arrays are automatically downloaded from <a href="https://zenodo.org/records/13764770">zenodo.org/records/13764770</a> when Humatch is first run.

$${\color{red}If \space you \space have \space issues \space with \space the \space auto \space downloads}$$ then the 3 weights (.h5) files and 24 germline likeness lookup arrays (.npy) files can be manually downloaded and saved in ```Humatch/Humatch/trained_models``` and ```Humatch/Humatch/germline_likeness_lookup_arrays``` respectively. Once these files are downloaded and saved in the right folders, please rerun the ```pip install .``` command to add these files to Humatch's package data.

## Web tool

We have a live [web tool](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/humatch/):

<https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/humatch/>

If you use this code or the web tool in your work, please [cite the paper](#citation).

## Humanness classification

Humatch can be used to obtain heavy, light, and paired predictions for a given VH/VL pair. As Humatch was trained on complete VH/VL sequences, only complete sequences should be used as input. An example notebook is provided and predictions can also be obtained from the command line e.g.

```
Humatch-classify
    -H EVQLVESGGG...VSS
    -L DIVMTQGALP...EIK
    -s
```

Output:

```
hv:     hv3
lv:     kv2
CNN_H:  0.000
CNN_L:  0.000
CNN_P:  0.000
```

Predictions for all heavy and light V-genes are returned by ommitting the ```-s``` summary flag, otherwise only the top-scoring human v-genes are returned. Classification can be run for individual heavy/light chains - in this instance, only the heavy/light CNN score will be returned.

For high throughput screening, Humatch can be run on a csv file of antibody sequences e.g.

```
Humatch-classify
    -i data/example.csv
    --vh_col heavy
    --vl_col light
```

Output:

<div align="center">

| VH | VL | hv1 | ... | hv7 | lv1 | ... | lv10 | kv1 | ... | kv7 | CNN_P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| QVQ...VSS | EIV...IK- | 0.999 | ... | 0.000 | 0.000 | ... | 0.000 | 0.000 | ... | 0.000 | 0.998 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

</div>

The output csv will contain Humatch's predictions alongside the aligned, padded VH/VL sequences of the columns specified. Output paths can be specified with the ```-o``` argument.

## Humanisation

Humatch is primarily designed to offer experimental-like humanisation in seconds. Like humanness classification, an example notebook is provided in addition to the command line interface e.g.

```
Humatch-humanise
    -H QVNLLQSGAA...VSA
    -L DTVLTQSPAL...EIK
    -v
```

Output:

```
Humanised sequences:
        EVKLVESGGG...VSS
        DTVLTQSPAL...EIK
        Edit:   24
        HV:     hv1
        LV:     kv3
        CNN_H:  0.958
        CNN_L:  0.999
        CNN_P:  0.981
```

Using the verbose ```-v``` flag will show you the default config parameters used by Humatch. Users may design their own config file and point to this instead using the ```--config``` argument if they wish to specify target genes or add/remove residues Humatch cannot mutate.

If humanising many sequences, Humatch can be run on a csv file of antibody sequences similarly to classification e.g.

```
Humatch-humanise
    -i data/example.csv
    --vh_col heavy
    --vl_col light
```

Output (the first example sequence is predicted to be human, so no edits are suggested):

<div align="center">

| Humatch_H | Humatch_L | Edit | HV | LV | CNN_H | CNN_L | CNN_P |
| --- | --- | --- | --- | --- | --- | --- | --- |
| QVQ...VSS | EIV...IK- | 0 | hv1 | kv3 | 0.999 | 1.0 | 0.998 |
| ... | ... | ... | ... | ... | ... | ... | ... |

</div>

## Sequence alignment

Users can run sequence alignment in isolation to determine how ANARCI has numbered a sequence. With this information, users can then specify IMGT positions to remain fixed during humanisation e.g. add ```CNN_fixed_imgt_positions_H: ["9 ", "81A", "120 "]``` to the config file (note - the spaces where insertion codes are not present are required)

```
Humatch-align
    -H QVQLVQSGAE...VSS
    -L EIVLTQSPVT...EIK
```

Output

```
1       Q       E
2       V       I
3       Q       V
3A      -       -
4       L       L
5       V       T
...     ...     ...
126     V       I
127     S       K
128     S       -
```

For small speed increases, users may also wish to pre-align many sequences to avoid repeating this step if wanting to trial many different humanisation configs e.g.

```
Humatch-align
    -i data/example.csv
    --vh_col heavy
    --vl_col light
```

This can be run with the ```--imgt_cols``` flag to return unique columns for each IMGT position, otherwise only two columns are returned - padded VH and VL. If csvs are pre-aligned (without the ```--imgt_cols``` flag), the alignment step can be avoided during classification and humanisation by including the ```--aligned``` flag.

## Citation

```
@article{Chinery2024,
  title = {Humatch - fast, gene-specific joint humanisation of antibody heavy and light chains},
  author = {Lewis Chinery, Jeliazko R Jeliazkov, and Charlotte M Deane},
  journal = {bioRxiv},
  year = {2024},
  doi = {10.1101/2024.09.16.613210}
}
```
