This repo contains code for extracting the molecules in
https://data.hpc.imperial.ac.uk/resolve/?doi=4618
into an AtomLite_ database.


.. _AtomLite: https://atomlite.readthedocs.io

Why?
=====

Because the original format of the published data is an out-of-date format,
namely an stk_ JSON dump.

.. _stk: https://stk.readthedocs.io

How?
====

The easiest thing to do is

.. code-block:: bash

  pip install cage-json-extractor

Now you can download the files

* ``cages.tar.gz`` - https://data.hpc.imperial.ac.uk/resolve/?doi=4618&file=3&access=
* ``cage_prediction.db`` - https://data.hpc.imperial.ac.uk/resolve/?doi=4618&file=2&access=

And run

.. code-block:: bash

  tar xf cages.tar.gz
  cage_json_extractor cages/amines2aldehydes3.json cage_prediction.db amines2aldehydes3.db

Now if we want to extract all the shape persistent 4+6 cages we can run

.. code-block:: bash

  extract_cages amines2aldehydes3.db FourPlusSix --output_directory extracted_cages

This will create a folder ``extracted_cages`` which holds a sub-folder for every
shape persistent 4+6 cage in ``amines2aldehydes3.db``. In the sub-folder you will
find the ``.mol`` file of the cage and its building blocks.

Enjoy! (and sorry I deprecated the ``.json`` files)

=)
