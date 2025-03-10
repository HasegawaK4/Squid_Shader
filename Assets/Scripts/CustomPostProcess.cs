using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;

[VolumeComponentMenu("Post-processing/Custom/BlurEffect")]
public class BlurEffect : CustomPostProcessVolumeComponent, IPostProcessComponent
{
    public ClampedFloatParameter blurSize = new ClampedFloatParameter(1.0f, 0.0f, 10.0f);

    private Material material;

    public override void Setup()
    {
        if (material == null)
        {
            Shader shader = Shader.Find("Hidden/Custom/BlurShader");
            material = new Material(shader);
        }
    }

    public override void Render(CommandBuffer cmd, HDCamera camera, RTHandle source, RTHandle destination)
    {
        material.SetFloat("_BlurSize", blurSize.value);
        HDUtils.DrawFullScreen(cmd, material, destination, source);
    }

    public bool IsActive() => blurSize.value > 0.0f;
}