#!/usr/bin/python3

from os.path import join
from uuid import uuid4
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

class StructureParser(object):
  def __init__(self, cif_path):
    parser = MMCIFParser(QUIET = True)
    self.structure = parser.get_structure(str(uuid4()), cif_path)
    self.mmcif_dict = MMCIF2Dict(cif_path)
  def get_entities(self,):
    ents = dict(
      id = self.mmcif_dict.get("_struct_asym.id", []), # chain id of a entity
      entity_id = self.mmcif_dict.get("_struct_asym.entity_id", []), # entity id
      detail = self.mmcif_dict.get("_struct_asym.details", []), # chain name
    )
    ents = pd.DataFrame(ents)
    entities = dict()
    for _, row in ents.iterrows():
      if row['entity_id'] not in entities: entities[row['entity_id']] = list()
      entities[row['entity_id']].append({'chain_id': row['id'], 'detail': row['detail']})
    return entities
  def get_connections(self, entity_id):
    entities = self.get_entities()
    assert entity_id in entities
    chain_ids = set([chain['chain_id'] for chain in entities[entity_id]])
    conn = dict(
      id = self.mmcif_dict.get("_struct_conn.id", []), # connection id
      residue1_chain_id = self.mmcif_dict.get("_struct_conn.ptnr1_label_asym_id", []), # chain id
      residue1_name = self.mmcif_dict.get("_struct_conn.ptnr1_label_comp_id", []), # residue1 name
      residue1_seq_id = self.mmcif_dict.get("_struct_conn.ptnr1_label_seq_id", []), # residue1 seq id in chain
      residue2_chain_id = self.mmcif_dict.get("_struct_conn.ptnr2_label_asym_id", []), # chain id
      residue2_name = self.mmcif_dict.get("_struct_conn.ptnr2_label_comp_id", []), # residue2 name
      residue2_seq_id = self.mmcif_dict.get("_struct_conn.ptnr2_label_seq_id", []), # residue2 seq id in chain
    )
    conn = pd.DataFrame(conn)
    connections = conn[conn['residue1_chain_id'].isin(chain_ids)|conn['residue2_chain_id'].isin(chain_ids)]
    return connections
  def get_chain(self, chain_id):
    structure_id = str(uuid4()) if structure_id is None else structure_id
  def get_structure(self,):
    return self.structure

if __name__ == "__main__":
  parser = StructureParser(join('tests', '1sfi.cif'))
  entities = parser.get_entities()
  print(entities)
  connections = parser.get_connections('1')
  print(connections)
  structure = parser.get_structure()
  for model_id, model in enumerate(structure):
    # 对于一个构象
    output_model = []
    for chain in model:
      # 对于一个chain
      print(f"model_id: {model_id}, chain_id: {chain.id}")
      for residue in chain:
        # 对于一个residue
        print(residue.get_resname(), end = ' ')
      print('')
