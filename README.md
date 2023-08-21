# SARS-CoV-2 antiviral prediction: REDIAL-2020

Predictor of several endpoints related to Sars-CoV-2. It provides predictions for Live Virus Infectivity, Viral Entry, Viral Replication, In Vitro Infectivity and Human Cell Toxicity using a combination of three models. Consensus results are obtained by averaging the prediction for the three different models for each activity and toxicity models. The models have been built using NCATS COVID19 data. Further details on result interpretations can be found here: https://drugcentral.org/Redial

## Identifiers

* EOS model ID: `eos8fth`
* Slug: `redial-2020`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Classification`
* Output: `Probability`
* Output Type: `Float`
* Output Shape: `Single`
* Interpretation: The model returns the probability of 1 (active) in each assay. Good drugs are active in CPE, 3CL and are inactive in cytotox, hCYTOX and ACE2 and/or are active in at least one of the following: AlphaLISA, CoV-PPE, MERS-PPE, while inactive in the counter screen, respectively: TruHit, CoV-PPE_cs, MERS-PPE_cs.

## References

* [Publication](https://www.nature.com/articles/s42256-021-00335-w#Sec9)
* [Source Code](https://github.com/sirimullalab/redial-2020/tree/v1.0)
* Ersilia contributor: [Pradnya2203](https://github.com/Pradnya2203)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos8fth)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8fth.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos8fth) (AMD64)

## Citation

If you use this model, please cite the [original authors](https://www.nature.com/articles/s42256-021-00335-w#Sec9) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a MIT license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!