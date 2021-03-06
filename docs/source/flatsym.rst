
======================================
Flatness/Symmetry module documentation
======================================

Overview
--------

The Flatness & Symmetry module (``flatsym``) allows a physicist to check their linac's beam profile
against well-known flatness/symmetry calculation standards, reporting absolute values. Film or EPID images
can be loaded in and analyzed.

Running the Demo
----------------

To run the demo of the ``flatsym`` module, import the main class from it and run the demo method::

    from pylinac.flatsym import BeamImage

    my_img = BeamImage()
    my_img.run_demo()

Which will result in the following plot:

.. image:: images/flatsym_demo.png

Typical Use
-----------

In most instances, a physicist is interested in quickly calculating the flatness, symmetry, or both of the
image in question. The ``flatsym`` module allows you to do this easily, using any of multiple definitions of flatness
or symmetry.

To get started, import the ``BeamImage`` class from the ``flatsym`` module::

    from pylinac.flatsym import BeamImage

Loading images is easy and just like any other module::

    # from a file
    my_file = r"C:/my/QA/folder/img.dcm"
    my_img = BeamImage(my_file)

    # or using a UI
    my_img = BeamImage.open_UI()

You can calculate the flatness or symmetry directly, e.g.::

    my_img.flatness()
    [1.91, 2.6]  # or whatever the values are. Returns: [crossplane, inplane]
    my_img.symmetry()
    [3.08, 2.3]

Any physicist worth their salt will want to check that the right profile has been taken however.
This is easily done with similar methods::

    my_img.plot_flatness()

.. image:: images/flat_demo.png

Default behavior is to analyze both planes at the calculated center of the field, but other options exist:

* Individual planes can be analyzed instead of both.
* The position of the profile can be directly passed.
* Multiple :ref:`analysis_definitions` exist.

.. seealso:: :meth:`~pylinac.flatsym.BeamImage.plot_flatsym()` parameters for more details.

.. _analysis_definitions:

Analysis Definitions
--------------------

.. warning:: The following definitions are for **photons** only.

There are multiple definitions for both flatness and symmetry. Your machine vendor uses certain equations,
or your clinic may use a specific definition. Pylinac has a number of built-in definitions which you can use.

Symmetry:

* -- Name, Vendors that use it -- Equation
* -- **Point Difference, Varian** -- :math:`100 * max(|L_{pt} - R_{pt}|)/ D_{CAX}` over 80%FW, where :math:`L_{pt}` and :math:`R_{pt}` are
  equidistant from the CAX.
* -- **Point Difference Quotient (IEC), Elekta** -- :math:`100 * max(|L_{pt}/R_{pt}|, |R_{pt}/L_{pt}|)` over 80%FW if 10<FW<30cm [#elekta]_.

Flatness:

* -- Name, Vendors that use it -- Equation
* -- **Variation over mean (80%), Varian** -- :math:`100 * |D_{max} - D_{min}| / (D_{max} + D_{min})` within 80%FW.
* -- **Dmax/Dmin (IEC), Elekta** -- :math:`100 * D_{max}/D_{min}` within 80%FW for 10<FW<30cm [#elekta]_.

.. note:: Siemens and other definitions (e.g. Area, Area/2) will be added if the community `asks <https://github.com/jrkerns/pylinac/issues>`_ for it.

.. [#elekta]
    The region calculated over actually varies by the following: for 5<FW<10cm, FW - 2*1cm; for 10<FW<30cm, FW - 2*0.1*FW (i.e. 80%FW);
    for 30cm<FW, FW - 2*6cm. Pylinac currently only uses the 80%FW no matter the FW, but accounting for the FW will come in a future
    version.

Algorithm
---------

There is little of a true "algorithm" in ``flatsym`` other than automatic field determination. Thus, this section is more terminology and
notekeeping.

**Allowances**

* The image can be any size.
* The image can be digitized film or EPID (most image formats and DICOM).
* The image can be either inversion (Radiation is dark or light).

**Restrictions**

* The module is only meant for photon analysis at the moment (there are sometimes different equations for electrons for the same
  definition name).
* The image should be near perpendicular/normal to the image edge; this actually won't cause the module to fail, but may invalidate the
  results' accuracy. Accounting for angle may come in a future version if the community asks for it.

**Analysis**

* *Determine the profile position* - If the profile position determination is left to the automatic analysis, the row and column image sums
  are analyzed
  for the center of the FWHM. These then determine the location, either crossplane, inplane, or both. If the values are explicitly passed
  in, these values are used directly (or converted from a fraction).
* *Extract profiles* - With the positions known, profile(s) are extracted and analyzed according to the method specified (see
  :ref:`analysis_definitions`). For symmetry calculations that operate around the CAX, the CAX must first be determined, which is
  the center of the FWHM of the profile.

API Documentation
-----------------

.. autoclass:: pylinac.flatsym.BeamImage
    :no-show-inheritance:


