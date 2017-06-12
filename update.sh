#!/bin/bash
cp /cloud/Code/pyapp/bp.py /cloud/Code/pycomplete/pybp
sed -i -- 's/include = os.environ/#include = os.environ/g' /cloud/Code/pycomplete/pybp/bp.py
sed -i -- 's/sys.path.append/#sys.path.append/g' /cloud/Code/pycomplete/pybp/bp.py
cp /cloud/Code/pylib/*.py /cloud/Code/pycomplete/pybp
