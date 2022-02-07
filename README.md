# YRBS Predictions

#### R code

| File (Markdown with output)                                        | File (R code)                                                        | Purpose                                                                                                                               |
|--------------------------------------------------------------------|----------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| [train_and_predict.md](train_and_predict.md)                       | [train_and_predict.Rmd](train_and_predict.Rmd)                       | Generates LOO (leave-one-out) predictions for each state for which we have observed responses to one or both of the focal questions.  |
| [train_and_predict_male_q66.md](train_and_predict_male_q66.md)     | [train_and_predict_male_q66.Rmd](train_and_predict_male_q66.Rmd)     | Same as `train_and_predict` but just for male respondents and the same-sex contact question.                                          |
| [evaluate_loo_predictions.md](evaluate_loo_predictions.md)         | [evaluate_loo_predictions.Rmd](evaluate_loo_predictions.Rmd)         | Evaluates the LOO predictions against the true proportions we observe in YRBS data.                                                   |
| [make_final_state_predictions.md](make_final_state_predictions.md) | [make_final_state_predictions.Rmd](make_final_state_predictions.Rmd) | Combines the LOO predictions with true proportions and makes predictions for all the states without responses to the focal questions. |

#### Other contents
* `data/` contains input and output CSVs used in the R code. NOTE: The raw survey data is not public data. Contact the author(s) to get instructions for how to acquire the raw datasets.

#### Requirements
* R packages loaded into the R notebooks
