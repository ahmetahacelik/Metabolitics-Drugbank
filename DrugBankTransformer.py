from sklearn_utils.utils import SkUtilsIO
from sklearn_utils.utils import feature_importance_report
from sklearn.decomposition import PCA
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
import ast
from sklearn.pipeline import Pipeline

class Drugbank:
    def __init__(self,file,file2):
        self.file=file
        self.file2=file2

    def pathway(self,file):
        self.X=[]
        self.all=[]
        for line in self.file:
            line = line.rstrip()
            self.all.append(line)
        for i in range(214):
            my_dict = ast.literal_eval(self.all[i])
            self.X.append(my_dict)
        return self.X

    def normal_pathway(self,file2):
        self.X2=[]
        self.all2=[]
        for line in self.file2:
            line=line.rstrip()
            self.all2.append(line)
        for i in range(214):
            my_dict = ast.literal_eval(self.all2[i])
            self.X2.append(my_dict)

        return self.X2
    def transformer(self,x,y):

        pipe = Pipeline([
            ('vect', DictVectorizer(sparse=False)),
            ('pca', PCA()),
            ('clf', LogisticRegression(C=0.3e-6, random_state=43, solver='lbfgs'))
        ])

        kf = StratifiedKFold(n_splits=10, random_state=43)

        scores = cross_val_score(pipe, x, y, cv=kf, n_jobs=-1, scoring='f1_micro')
        print('kfold test: %s' % scores)
        print('mean: %s' % scores.mean().round(3))
        print('std: %s' % scores.std().round(3))

    def feature_importance_report(self,x,y):

        df = feature_importance_report(x, y)
        print(df)


if __name__ == "__main__":
    x, y = SkUtilsIO('BC.csv').from_csv(label_column='stage')
    y = ['healthy' if i == 'h' else 'bc' for i in y]
    file = open("pathwaying", "r")
    file2 = open("normal_pathway", "r")
    bank=Drugbank(file,file2)
    bank.pathway(file)
    bank.transformer(bank.normal_pathway(file2),y)
