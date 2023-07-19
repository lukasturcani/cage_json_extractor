import argparse
import pathlib
from typing import cast

import atomlite
import rdkit.Chem.AllChem as rdkit


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(exist_ok=True, parents=True)
    db = atomlite.Database(args.database)
    num_collapsed = 0
    num_persistent = 0
    num_undetermined = 0
    for entry in db.get_entries():
        if (
            "inchi_building_blocks" in entry.properties
            and entry.properties["collapsed"] is False
            and args.topology in cast(str, entry.properties["topology"])
        ):
            num_persistent += 1
            smiles_building_blocks = cast(
                list[str],
                entry.properties["smiles_building_blocks"],
            )
            inchi_building_blocks = cast(
                list[str],
                entry.properties["inchi_building_blocks"],
            )

            cage_name = "_".join(smiles_building_blocks)
            cage_directory = args.output_directory / cage_name
            cage_directory.mkdir(exist_ok=True, parents=True)

            cage = atomlite.json_to_rdkit(entry.molecule)
            building_blocks = [
                atomlite.json_to_rdkit(bb_entry.molecule)
                for bb_entry in db.get_entries(inchi_building_blocks)
            ]

            rdkit.MolToMolFile(
                cage,
                str(cage_directory / "cage.mol"),
                forceV3000=True,
            )
            for i, building_block in enumerate(building_blocks):
                rdkit.MolToMolFile(
                    building_block,
                    str(cage_directory / f"bb_{i}.mol"),
                    forceV3000=True,
                )
        elif (
            "inchi_building_blocks" in entry.properties
            and entry.properties["collapsed"] is True
            and args.topology in cast(str, entry.properties["topology"])
        ):
            num_collapsed += 1
        elif (
            "inchi_building_blocks" in entry.properties
            and entry.properties["collapsed"] is None
            and args.topology in cast(str, entry.properties["topology"])
        ):
            num_undetermined += 1

    print(f"Total cages: {num_collapsed+num_persistent+num_undetermined}")
    print(f"Num collapsed: {num_collapsed}")
    print(f"Num persistent: {num_persistent}")
    print(f"Num undetermined: {num_undetermined}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=pathlib.Path)
    parser.add_argument("topology", type=pathlib.Path)
    parser.add_argument(
        "--output_directory",
        type=pathlib.Path,
        default=pathlib.Path("cages"),
    )
    return parser.parse_args()
