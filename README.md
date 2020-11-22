# Machine Learning - Drug Resistance Classification 
The presence/absence of 4,016 genes was determined for 85 *M. tuberculosis* clinical isolates, all of which had been previously labeled with a drug resistance (DR) classification of either mono-resistant, multidrug resistant, or extensively drug resistant.

To create models that learn to associate features (genes) with labels (DR classifications), 68/85 isolates (80%) were used to train the various machine learning models. 17/85 isolates (20%) were then used as a test set for which model accuracy was calculated. All machine learning model accuracies were compared to one another to evaluate their performance side-by-side.

## Summary
Machine Learning Model | Accuracy | Rank
--- | --- | ---
MLP Neural Network (3 layers + 2 dropouts) | 76% | 1st-t
Random Forest (10,000 trees) | 76% | 1st-t
SVM (RBF kernel) | 76% | 1st-t
Decision Tree | 71% | 4th-t
SVM (poly kernel) | 71% | 4th-t
MLP Neural Network (3 layers) | 65% | 6th-t
Logistic Regression | 65% | 6th-t

The MLP neural network with dropout layers, random forest classifier, and SVM classifier (RBF kernel) tied for first with accuracies of 76%. This means that each model was able to correctly predict DR classifications for 13/17 *M. tuberculosis* isolates in the test set based on the specific gene presence/absence profile for each isolate. With more data (genomes) and more parameter testing, it's likely a clear-cut winner would emerge.

## Requirements 
Python 3.7.4

A complete list of Python dependencies can be found in ```requirements.txt```.