import argparse
import json
import pathlib
import sqlite3

import atomlite
import rdkit.Chem.AllChem as rdkit
import stk


def main() -> None:
    args = parse_args()

    db = atomlite.Database(args.output_database)
    property_db = sqlite3.connect(args.property_database)

    with open(args.json) as f:
        json_db = json.load(f)

    for cage_json in json_db:
        if "FourPlusSix" not in cage_json["topology"]:
            continue
        smiles_building_blocks: list[atomlite.Json] = []
        inchi_building_blocks: list[atomlite.Json] = []
        for bb_json in cage_json["building_blocks"]:
            ((_, bb_mol_block),) = bb_json["conformers"]
            bb = rdkit.MolFromMolBlock(
                bb_mol_block,
                sanitize=False,
                removeHs=False,
            )
            stk_bb = stk.BuildingBlock.init_from_rdkit_mol(bb)
            inchi = stk.Inchi().get_key(stk_bb)
            smiles = stk.Smiles().get_key(stk_bb)
            db.update_entries(
                entries=atomlite.Entry.from_rdkit(
                    key=inchi,
                    molecule=bb,
                    properties={
                        "smiles": smiles,
                        "inchi": inchi,
                    },
                ),
                commit=False,
            )
            smiles_building_blocks.append(smiles)
            inchi_building_blocks.append(inchi)

        ((_, cage_mol_block),) = cage_json["conformers"]
        cage = rdkit.MolFromMolBlock(
            cage_mol_block,
            sanitize=False,
            removeHs=False,
        )
        stk_cage = stk.BuildingBlock.init_from_rdkit_mol(cage)
        db.update_entries(
            entries=atomlite.Entry.from_rdkit(
                key=stk.Inchi().get_key(stk_cage),
                molecule=cage,
                properties={
                    "smiles_building_blocks": smiles_building_blocks,
                    "inchi_building_blocks": inchi_building_blocks,
                    "name": cage_json["name"],
                    "collapsed": get_collapsed(property_db, cage_json["name"]),
                },
            ),
            commit=False,
        )

    db.connection.commit()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("json", type=pathlib.Path)
    parser.add_argument("property_database", type=pathlib.Path)
    parser.add_argument("output_database", type=pathlib.Path)
    return parser.parse_args()


def get_collapsed(db: sqlite3.Connection, name: str) -> bool | None:
    (value,) = db.execute(
        "SELECT collapsed "
        "FROM cages "
        "WHERE name=? "
        "AND reaction LIKE '%amine2aldehyde3%' "
        "AND topology LIKE '%FourPlusSix%'",
        (name,),
    ).fetchone()
    if value == 0 or value == 1:
        return bool(value)
    elif value is None:
        return None
    else:
        raise RuntimeError("goofed")