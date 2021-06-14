.ONESHELL:
SHELL := /bin/bash
SRC = $(wildcard ./*.ipynb)

all: pyndamics3 docs

pyndamics3: $(SRC)
	nbdev_build_lib
	touch pyndamics3

sync:
	nbdev_update_lib

docs_serve: docs
	cd docs && bundle exec jekyll serve

docs: $(SRC)
	nbdev_build_docs
	touch docs

github: pyndamics3
	open -a Github\ Desktop

install: pyndamics3
	pip install . --upgrade

test:
	nbdev_test_nbs

release: pypi
	nbdev_conda_package
	nbdev_bump_version

pypi: dist
	twine upload --repository pypi dist/*

dist: clean
	python setup.py sdist bdist_wheel

clean:
	rm -rf dist