# CTHMMs-for-adaptive-sampling

Repository for a study proposing a method for adaptive sampling in Ecological Momentary Assessment type settings, using continuous-time hidden Markov models. 

Part of the EPSRC-funded [Wearable Clinic](https://www.herc.ac.uk/research_project/wearable-clinic-connecting-health-self-care/) project. 

## Abstract

Wearable and mobile technology provides new opportunities to manage health conditions remotely and unobtrusively.
For example, healthcare providers can repeatedly sample a person's condition to monitor disease progression and intervene if necessary.
To ensure that this technology is deployed effectively, there is a utility-tolerability trade-off between collecting information at sufficient frequencies and quantities to be useful, and over-burdening the user or the underlying technology. 
Selecting the next sampling time adaptively using previous responses, so that people are only sampled at high frequency when necessary, can help to manage this trade-off.

We present an approach to adaptive sampling using clusters continuous-time hidden Markov models.
The model is used to predict, at any given sampling time, the probability of moving to a high-risk state, and the next sample time is scheduled when this probability has exceeded a given threshold.
The clusters, each representing a distinct sub-model, allow heterogeneity in states and state transitions. 

The work is illustrated using longitudinal mental-health symptom data in 49 people collected using ClinTouch, a mobile app designed to monitor patients with schizophrenia. 
Using these data we show how the adaptive sampling scheme behaves under different model parameters and risk thresholds, and how, subject to accepting a slightly higher average detection delay, the sampling frequency can be substantially reduced on average. 