[&larr; Contents](./__contents.md)

---

# Distributing Package

Do not forget update version at ``odelabs/__init__.py``

1. Create package:
    ```shell script
    # python setup.py sdist -d ./__dist
    python setup.py egg_info --egg-base ./__dist sdist -d ./__dist
    ```

1. Upload it to pypi (you must provide your pypi login and account):
    ```shell script
    twine upload __dist/odelabs-A.B.C
    ```

If you have any doubts, visit:

* [Versioning Scheme](https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/specification.html)
* [Packaging and distributing projects](https://packaging.python.org/guides/distributing-packages-using-setuptools/)
* [Writing the Setup Script](https://docs.python.org/3/distutils/setupscript.html)
* [pypi Classifiers](https://pypi.org/classifiers/)
* [How to upload your python package to pypi](https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56)

---

[&larr; Contents](./__contents.md)