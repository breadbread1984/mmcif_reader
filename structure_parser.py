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
  def get_residues(self,):
    poly = dict(
      chain_id = self.mmcif_dict.get("_pdbx_poly_seq_scheme.asym_id", []), # chain_id
      seq_id = self.mmcif_dict.get("_pdbx_poly_seq_scheme.seq_id", []), # seq id in chain
      residue_name = self.mmcif_dict.get("_pdbx_poly_seq_scheme.mon_id", []), # residue name
    )
    poly = pd.DataFrame(poly)
    poly['seq_id'] = poly['seq_id'].astype(int)
    nonpoly = dict(
      chain_id = self.mmcif_dict.get("_pdbx_nonpoly_scheme.asym_id", []), # chain_id
      suedoresidue_name = self.mmcif_dict.get("_pdbx_nonpoly_scheme.mon_id", []), # pseudo residue name
    )
    nonpoly = pd.DataFrame(nonpoly)
    return poly, nonpoly
  def get_structure(self,):
    return self.structure

if __name__ == "__main__":
  parser = StructureParser(join('tests', '1sfi.cif'))
  entities = parser.get_entities()
  print(entities)
  connections = parser.get_connections('1')
  print(connections)
  poly, nonpoly = parser.get_residues()
  # print chain A
  print(poly[poly['chain_id'] == 'A'].sort_values('seq_id', ascending = True))