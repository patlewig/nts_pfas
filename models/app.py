from model_app import MoleculeProcessor, CategoryManager, CategoryPrediction,categories

category_manager = CategoryManager(categories)

rev_dict = 'reverse_dict2.pkl'
rf_model ='final_model_v2.sav'

df = category_manager.make_df('CCCN(CCNC(=O)c1ccc(Cc2ccc(C(O)=O)cc2)cc1)S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', dtx = '123', cat="Aromatic PFASs")

catpred = CategoryPrediction(rev_dict, rf_model)

print(catpred.make_prediction(df))