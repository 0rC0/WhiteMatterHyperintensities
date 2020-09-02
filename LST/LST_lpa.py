"""
Nipype module for Lesion segmentation toolbox lesion prediction algorithm (LST_lpa)
"""

from nipype.interfaces.base import TraitedSpec, \
    BaseInterface, BaseInterfaceInputSpec, File, traits
from nipype.interfaces.matlab import MatlabCommand, MatlabInputSpec
import os
from string import Template


class LPAInputSpec(MatlabInputSpec):
    flair = File(exists=True, mandatory=True)
    #anat = File(exists=True, mandatory=False)
    #report = traits.Enum('0', '1', usedefault=True, mandatory=False)


class LPAOutputSpec(TraitedSpec):
    matlab_output = traits.Str()
    #out_file = File()


class LPA(MatlabCommand):
    """
    Here is needed a description

    ToDo:
        implement T1w and bool Report as input
        implemet output. Should contain:
            -) lesion mask
            -) mat file
            -) HTML report

    Example:
    >>> import LST_lpa
    >>> lpa=LST_lpa.LPA()
    >>> lpa.inputs.flair='/home/orco/WMH/WMHSP_test20200830/sub-SMART008_FLAIR.nii.gz'
    >>> out=lpa.run()
    """
    input_spec = LPAInputSpec
    output_spec = LPAOutputSpec

    def _my_script(self):
        script = """
        ps_LST_lpa('%s','',0)
        """ % self.inputs.flair
        return script

    def run(self, **inputs):
        # Inject your script
        self.inputs.script = self._my_script()
        results = super(MatlabCommand, self).run(**inputs)
        stdout = results.runtime.stdout
        # Attach stdout to outputs to access matlab results
        results.outputs.matlab_output = stdout
        return results

    def _list_outputs(self):
        outputs = self._outputs().get()
        return outputs