# SARS-CoV-2 antiviral prediction: REDIAL-2020

Predictor of several endpoints related to Sars-CoV-2. It provides predictions for Live Virus Infectivity, Viral Entry, Viral Replication, In Vitro Infectivity and Human Cell Toxicity using a combination of three models. Consensus results are obtained by averaging the prediction for the three different models for each activity and toxicity models. The models have been built using NCATS COVID19 data. Further details on result interpretations can be found here: https://drugcentral.org/Redial

This model was incorporated on 2023-03-27.

## Information
### Identifiers
- **Ersilia Identifier:** `eos8fth`
- **Slug:** `redial-2020`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `COVID-19`
- **Target Organism:** `SARS-CoV-2`
- **Tags:** `Sars-CoV-2`, `COVID19`, `Antiviral activity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `11`
- **Output Consistency:** `Fixed`
- **Interpretation:** The model returns the probability of the molecule being active in each assay. Good drugs: active in CPE, 3CL and inactive in cytotox, hCYTOX and ACE2 and/or active in at least one of: AlphaLISA, CoV-PPE, MERS-PPE while inactive in the counter screen respectively: TruHit, CoV-PPE_cs, MERS-PPE_cs.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| 3cl | float | high | Probability that the compound inhibits SARS-CoV-2 3CL (Mpro) |
| ace2 | float | high | Probability that the compound inhibits ACE2 directly causing secondary effects |
| alphalisa | float | high | Probability that a compound disrupts the Spike-ACE2 protein-protein interaction |
| cpe | float | high | CytoPhatic Effect or probability that a compound reverses the virus effect on Vero E6 cells |
| cov1_ppe | float | high | Pprobability that the compound inhibits viral entry into cells using Sars-Cov1 pseudotype particles |
| cov1_ppe_cs | float | high | Counterscreen for COV1-PPE |
| mers_ppe | float | high | Counterscreen for MERS-PPE |
| mers_ppe_cs | float | high | Probability that the compound inhibits viral entry into cells using MERS pseudotype particles |
| truhit | float | high | Probability that the AlpaLISA result is a false positive |
| cytotox | float | high | Probability of toxicity on Vero E6 cells as counterscreen for the CPE |

_10 of 11 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos8fth](https://hub.docker.com/r/ersiliaos/eos8fth)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8fth.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8fth.zip)

### Resource Consumption
- **Model Size (Mb):** `244`
- **Environment Size (Mb):** `1086`


### References
- **Source Code**: [https://github.com/sirimullalab/redial-2020/tree/v1.0](https://github.com/sirimullalab/redial-2020/tree/v1.0)
- **Publication**: [https://www.nature.com/articles/s42256-021-00335-w#Sec9](https://www.nature.com/articles/s42256-021-00335-w#Sec9)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2021`
- **Ersilia Contributor:** [Abellegese](https://github.com/Abellegese)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos8fth
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos8fth
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
