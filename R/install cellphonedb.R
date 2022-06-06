conda create --name cellphonedb python=3.5.5
conda activate cellphonedb
pip install cellphonedb

vi /home/yah2014/anaconda3/envs/cellphonedb/lib/python3.5/site-packages/cellphonedb/src/core/database/sqlalchemy_repository/InteractionRepository.py

change
"multidata_expanded: pd.DataFrame = self.database_manager.get_repository('multidata').get_all_expanded(include_gene)"
to
"multidata_expanded = self.database_manager.get_repository('multidata').get_all_expanded(include_gene)"
