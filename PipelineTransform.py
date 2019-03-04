from sklearn.preprocessing import LabelEncoder
from metabolitics.preprocessing import MetaboliticsPipeline
from sklearn_utils.utils import SkUtilsIO
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
le=LabelEncoder()
X, y = SkUtilsIO('BC.csv').from_csv(
    label_column='stage')
y = ['healthy' if i == 'h' else  'bc' for i in y]

pre=  MetaboliticsPipeline([
        'metabolite-name-mapping',
        'standard-scaler',
        'metabolitics-transformer',
        'reaction-diff',
        'feature-selection',
        'pathway-transformer',
        'transport-pathway-elimination'
    ])
X_fs_pathways = pre.fit_transform(X, y)


